#include "TruthTrackerAlg.h"
#include "DataHelper/HelixClass.h"
#include "DDSegmentation/Segmentation.h"
#include "DetSegmentation/GridDriftChamber.h"
#include "DetInterface/IGeomSvc.h"
#include "DD4hep/Detector.h"
#include "DD4hep/IDDescriptor.h"
#include "DD4hep/Plugins.h"
#include "DD4hep/DD4hepUnits.h"
#include "UTIL/ILDConf.h"

//external
#include "CLHEP/Random/RandGauss.h"
#include "edm4hep/EventHeaderCollection.h"
#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/SimTrackerHitCollection.h"
#include "edm4hep/TrackerHitCollection.h"
#include "edm4hep/TrackCollection.h"
#include "edm4hep/MCRecoTrackerAssociationCollection.h"
#include "edm4hep/MCRecoParticleAssociationCollection.h"
#include "edm4hep/ReconstructedParticleCollection.h"

DECLARE_COMPONENT(TruthTrackerAlg)

TruthTrackerAlg::TruthTrackerAlg(const std::string& name, ISvcLocator* svcLoc)
: GaudiAlgorithm(name, svcLoc),m_dd4hep(nullptr),m_gridDriftChamber(nullptr),
    m_decoder(nullptr)
{
    declareProperty("NoiseSimHitsCollection", w_SimNoiseHCol,
            "Handle of the output simulation  DC Noise hits collection");
    declareProperty("NoiseDCHitsCollection", w_NoiseHitCol,
            "Handle of the input DC Noise hits collection");
    declareProperty("NoiseDCHitAssociationCollection", w_NoiseAssociationCol,
            "Handle of the simNoiseDCHit and noiseDC Noise hits collection");
    declareProperty("MCParticle", m_mcParticleCol,
            "Handle of the input MCParticle collection");
    declareProperty("DriftChamberHitsCollection", m_DCSimTrackerHitCol,
            "Handle of DC SimTrackerHit collection");
    declareProperty("DigiDCHitCollection", m_DCDigiCol,
            "Handle of DC digi(TrackerHit) collection");
    declareProperty("DCHitAssociationCollection", m_DCHitAssociationCol,
            "Handle of association collection");
    declareProperty("DCTrackCollection", m_DCTrackCol,
            "Handle of DC track collection");
    declareProperty("SubsetTracks", m_siSubsetTrackCol,
            "Handle of silicon subset track collection");
    declareProperty("SDTTrackCollection", m_SDTTrackCol,
            "Handle of SDT track collection");
    declareProperty("VXDTrackerHits", m_VXDTrackerHits,
            "Handle of input VXD tracker hit collection");
    declareProperty("SITTrackerHits", m_SITTrackerHits,
            "Handle of input SIT tracker hit collection");
    declareProperty("SETTrackerHits", m_SETTrackerHits,
            "Handle of input SET tracker hit collection");
    declareProperty("FTDTrackerHits", m_FTDTrackerHits,
            "Handle of input FTD tracker hit collection");
    declareProperty("SITSpacePoints", m_SITSpacePointCol,
            "Handle of input SIT hit collection");
//    declareProperty("SETSpacePoints", m_SETSpacePointCol,
//            "Handle of input SET hit collection");
    declareProperty("FTDSpacePoints", m_FTDSpacePointCol,
            "Handle of input FTD hit collection");
    declareProperty("VXDCollection", m_VXDCollection,
            "Handle of input VXD hit collection");
    declareProperty("SITCollection", m_SITCollection,
            "Handle of input SIT hit collection");
    declareProperty("SETCollection", m_SETCollection,
            "Handle of input SET hit collection");
    declareProperty("FTDCollection", m_FTDCollection,
            "Handle of input FTD hit collection");
    declareProperty("TruthTrackerHitCollection", m_truthTrackerHitCol,
            "Handle of output truth TrackerHit collection");
    declareProperty("SmearTrackerHitCollection", m_SmeartruthTrackerHitCol,
            "Handle of output smear truth TrackerHit collection");
    declareProperty("SmearSimHitsCollection", w_SimSmearHCol,
            "Handle of output smear SimTrackerHit collection");
    declareProperty("SmearDCHitAssociationCollection", w_SmearAssociationCol,
            "Handle of output smear simulationhit and TrackerHit collection");
}


StatusCode TruthTrackerAlg::initialize()
{
    ///Get geometry
    m_geomSvc=service<IGeomSvc>("GeomSvc");
    if (!m_geomSvc) {
        error()<<"Failed to get GeomSvc."<<endmsg;
        return StatusCode::FAILURE;
    }
    ///Get Detector
    m_dd4hep=m_geomSvc->lcdd();
    if (nullptr==m_dd4hep) {
        error()<<"Failed to get dd4hep::Detector."<<endmsg;
        return StatusCode::FAILURE;
    }
    //Get Field
    m_dd4hepField=m_geomSvc->lcdd()->field();
    ///Get Readout
    dd4hep::Readout readout=m_dd4hep->readout(m_readout_name);
    ///Get Segmentation
    m_gridDriftChamber=dynamic_cast<dd4hep::DDSegmentation::GridDriftChamber*>
        (readout.segmentation().segmentation());
    if(nullptr==m_gridDriftChamber){
        error()<<"Failed to get the GridDriftChamber"<<endmsg;
        return StatusCode::FAILURE;
    }
    ///Get Decoder
    m_decoder = m_geomSvc->getDecoder(m_readout_name);
    if (nullptr==m_decoder) {
        error()<<"Failed to get the decoder"<<endmsg;
        return StatusCode::FAILURE;
    }

    m_tuple=nullptr;
    ///book tuple
    if(m_hist){
        NTuplePtr nt(ntupleSvc(), "TruthTrackerAlg/truthTrackerAlg");
        if(nt){
            m_tuple=nt;
        }else{
            m_tuple=ntupleSvc()->book("TruthTrackerAlg/truthTrackerAlg",
                    CLID_ColumnWiseTuple,"TruthTrackerAlg");
            if(m_tuple){
                StatusCode sc;
                sc=m_tuple->addItem("run",m_run);
                sc=m_tuple->addItem("evt",m_evt);
                sc=m_tuple->addItem("siMom",3,m_siMom);
                sc=m_tuple->addItem("siPos",3,m_siPos);
                sc=m_tuple->addItem("mcMom",3,m_mcMom);

                sc=m_tuple->addItem("mcPos",3,m_mcPos);
                sc=m_tuple->addItem("nDCTrackHit",m_nDCTrackHit,0,100000);
                sc=m_tuple->addItem("DriftDistance",m_nDCTrackHit,m_DriftDistance);
                sc=m_tuple->addItem("nSmearDCTrackHit",m_nSmearDCTrackHit,0,100000);
                sc=m_tuple->addItem("SmearDriftDistance",m_nSmearDCTrackHit,m_SmearDriftDistance);
                sc=m_tuple->addItem("nNoiseDCTrackHit",m_nNoiseDCTrackHit,0,100000);
                sc=m_tuple->addItem("NoiseDriftDistance",m_nNoiseDCTrackHit,m_NoiseDriftDistance);

                sc=m_tuple->addItem("nSimTrackerHitVXD",m_nSimTrackerHitVXD);
                sc=m_tuple->addItem("nSimTrackerHitSIT",m_nSimTrackerHitSIT);
                sc=m_tuple->addItem("nSimTrackerHitSET",m_nSimTrackerHitSET);
                sc=m_tuple->addItem("nSimTrackerHitFTD",m_nSimTrackerHitFTD);
                sc=m_tuple->addItem("nSimTrackerHitDC",m_nSimTrackerHitDC);
                sc=m_tuple->addItem("nTrackerHitVXD",m_nTrackerHitVXD);
                sc=m_tuple->addItem("nTrackerHitSIT",m_nTrackerHitSIT);
                sc=m_tuple->addItem("nTrackerHitSET",m_nTrackerHitSET);
                sc=m_tuple->addItem("nTrackerHitFTD",m_nTrackerHitFTD);
                sc=m_tuple->addItem("nTrackerHitDC",m_nTrackerHitDC);
                sc=m_tuple->addItem("nSpacePointSIT",m_nSpacePointSIT);
                sc=m_tuple->addItem("nSpacePointSET",m_nSpacePointSET);
                sc=m_tuple->addItem("nSpacePointFTD",m_nSpacePointFTD);
                sc=m_tuple->addItem("nHitOnSiTkXVD",m_nHitOnSiTkVXD);
                sc=m_tuple->addItem("nHitOnSiTkSIT",m_nHitOnSiTkSIT);
                sc=m_tuple->addItem("nHitOnSiTkSET",m_nHitOnSiTkSET);
                sc=m_tuple->addItem("nHitOnSiTkFTD",m_nHitOnSiTkFTD);
                sc=m_tuple->addItem("nHitOnSdtTkVXD",m_nHitOnSdtTkVXD);
                sc=m_tuple->addItem("nHitOnSdtTkSIT",m_nHitOnSdtTkSIT);
                sc=m_tuple->addItem("nHitOnSdtTkSET",m_nHitOnSdtTkSET);
                sc=m_tuple->addItem("nHitOnSdtTkFTD",m_nHitOnSdtTkFTD);
                sc=m_tuple->addItem("nHitOnSdtTkDC",m_nHitOnSdtTkDC);
                sc=m_tuple->addItem("nHitSdt",m_nHitOnSdtTk);
                sc=m_tuple->addItem("nNoiseOnSdt",m_nNoiseOnSdtTk);
            }
        }
    }
    return GaudiAlgorithm::initialize();
}

StatusCode TruthTrackerAlg::execute()
{
    info()<<"In execute()"<<endmsg;

    ///Output DC Track collection
    edm4hep::TrackCollection* dcTrackCol=m_DCTrackCol.createAndPut();

    ///Output SDT Track collection
    edm4hep::TrackCollection* sdtTkCol=m_SDTTrackCol.createAndPut();

    ///Output Hit collection
    auto truthTrackerHitCol=m_truthTrackerHitCol.createAndPut();

    ///Retrieve MC particle(s)
    const edm4hep::MCParticleCollection* mcParticleCol=nullptr;
    mcParticleCol=m_mcParticleCol.get();
    //if(m_mcParticleCol.exist()){mcParticleCol=m_mcParticleCol.get();}
    if(nullptr==mcParticleCol){
        debug()<<"MCParticleCollection not found"<<endmsg;
        return StatusCode::SUCCESS;
    }
    ///Retrieve DC digi
    const edm4hep::TrackerHitCollection* digiDCHitsCol=nullptr;
    const edm4hep::SimTrackerHitCollection* dcSimHitCol
        =m_DCSimTrackerHitCol.get();
    const edm4hep::MCRecoTrackerAssociationCollection* assoHits
        = m_DCHitAssociationCol.get();
    if(m_useDC){
        if(nullptr==dcSimHitCol){
            debug()<<"DC SimTrackerHitCollection not found"<<endmsg;
        }else{
            debug()<<"DriftChamberHitsCollection size "
                <<dcSimHitCol->size()<<endmsg;
            if(m_tuple)m_nSimTrackerHitDC=dcSimHitCol->size();
        }
        digiDCHitsCol=m_DCDigiCol.get();
        m_nDCTrackHit = digiDCHitsCol->size();
        if(nullptr==digiDCHitsCol){
            debug()<<"TrackerHitCollection not found"<<endmsg;
        }else{
            debug()<<"digiDCHitsCol size "<<digiDCHitsCol->size()<<endmsg;
            if((int) digiDCHitsCol->size()>m_maxDCDigiCut){
                debug()<<"Track cut by m_maxDCDigiCut "<<m_maxDCDigiCut<<endmsg;
                return StatusCode::SUCCESS;
            }
        }
    }

    // make noise
    edm4hep::SimTrackerHitCollection* SimVec = w_SimNoiseHCol.createAndPut();
    edm4hep::MCRecoTrackerAssociationCollection* AssoVec = 
        w_NoiseAssociationCol.createAndPut();
    edm4hep::TrackerHitCollection* Vec = w_NoiseHitCol.createAndPut();

    ////TODO
    //Output MCRecoTrackerAssociationCollection collection
    //const edm4hep::MCRecoTrackerAssociationCollection*
    //    mcRecoTrackerAssociationCol=nullptr;
    //if(nullptr==mcRecoTrackerAssociationCol){
    //    log<<MSG::DEBUG<<"MCRecoTrackerAssociationCollection not found"
    //        <<endmsg;
    //    return StatusCode::SUCCESS;
    //}
    //mcRecoTrackerAssociationCol=m_mcRecoParticleAssociation.get();

    ///New SDT track
    edm4hep::MutableTrack sdtTk=sdtTkCol->create();

    int nVXDHit=0;
    int nSITHit=0;
    int nSETHit=0;
    int nFTDHit=0;
    int nDCHitDCTk=0;
    int nDCHitSDTTk=0;

    ///Create track with mcParticle
    edm4hep::TrackState trackStateMc;
    getTrackStateFromMcParticle(mcParticleCol,trackStateMc);
    if(m_useTruthTrack.value()||!m_useSi){sdtTk.addToTrackStates(trackStateMc);}

    if(m_useSi){
        ///Retrieve silicon Track
        const edm4hep::TrackCollection* siTrackCol=nullptr;
        //if(m_siSubsetTrackCol.exist()){
        siTrackCol=m_siSubsetTrackCol.get();
        if(nullptr==siTrackCol){
            debug()<<"SDTTrackCollection is empty"<<endmsg;
            if(!m_useTruthTrack.value()) return StatusCode::SUCCESS;
        }else{
            debug()<<"SiSubsetTrackCol size "<<siTrackCol->size()<<endmsg;
            if(!m_useTruthTrack.value()&&0==siTrackCol->size()){
                return StatusCode::SUCCESS;
            }
        }
        //}

        for(auto siTk:*siTrackCol){
            if(!m_useTruthTrack.value()){
                debug()<<"siTk: "<<siTk<<endmsg;
                edm4hep::TrackState siTrackStat=siTk.getTrackStates(0);//FIXME?
                sdtTk.addToTrackStates(siTrackStat);
                sdtTk.setType(siTk.getType());
                sdtTk.setChi2(siTk.getChi2());
                sdtTk.setNdf(siTk.getNdf());
                sdtTk.setDEdx(siTk.getDEdx());
                sdtTk.setDEdxError(siTk.getDEdxError());
                sdtTk.setRadiusOfInnermostHit(
                        siTk.getRadiusOfInnermostHit());
            }
            if(!m_useSiTruthHit){
                debug()<<"use Si hit on track"<<endmsg;
                nVXDHit=addHotsToTk(siTk,sdtTk,lcio::ILDDetID::VXD,"VXD",nVXDHit);
                nSITHit=addHotsToTk(siTk,sdtTk,lcio::ILDDetID::SIT,"SIT",nSITHit);
                nSETHit=addHotsToTk(siTk,sdtTk,lcio::ILDDetID::SET,"SET",nSETHit);
                nFTDHit=addHotsToTk(siTk,sdtTk,lcio::ILDDetID::FTD,"FTD",nFTDHit);
            }//end of loop over hits on siTk
        }//end of loop over siTk

        if(m_useSiTruthHit){
            ///Add silicon SimTrackerHit
            debug()<<"Add silicon SimTrackerHit"<<endmsg;
            nVXDHit=addSimHitsToTk(m_VXDCollection,truthTrackerHitCol,sdtTk,"VXD",nVXDHit);
            nSITHit=addSimHitsToTk(m_SITCollection,truthTrackerHitCol,sdtTk,"SIT",nSITHit);
            nSETHit=addSimHitsToTk(m_SETCollection,truthTrackerHitCol,sdtTk,"SET",nSETHit);
            nFTDHit=addSimHitsToTk(m_FTDCollection,truthTrackerHitCol,sdtTk,"FTD",nFTDHit);
        }else{
            ///Add reconstructed hit or digi
            debug()<<"Add VXD TrackerHit"<<endmsg;
            nVXDHit=addHitsToTk(m_VXDTrackerHits,sdtTk,"VXD digi",nVXDHit);
            nSITHit=addHitsToTk(m_SITTrackerHits,sdtTk,"SIT digi",nSITHit);
            if(m_useSiSpacePoint.value()){
                ///Add silicon SpacePoint
                debug()<<"Add silicon SpacePoint"<<endmsg;
                if(m_useSiSpacePoint){
                    nSITHit=addHitsToTk(m_SITSpacePointCol,sdtTk,"SIT sp",nSITHit);
                }
               // nSETHit=addHitsToTk(m_SETSpacePointCol,sdtTk,"SET sp",nSETHit);
                nFTDHit=addHitsToTk(m_FTDSpacePointCol,sdtTk,"FTD sp",nFTDHit);
            }else{
                ///Add silicon TrackerHit
                debug()<<"Add silicon TrackerHit"<<endmsg;
                nSITHit=addHitsToTk(m_SITTrackerHits,sdtTk,"SIT digi",nSITHit);
                nSETHit=addHitsToTk(m_SETTrackerHits,sdtTk,"SET digi",nSETHit);
                nFTDHit=addHitsToTk(m_FTDTrackerHits,sdtTk,"FTD digi",nFTDHit);
            }//end of use space point
        }
    }//end of use silicon

    if(m_useDC){
        ///Create DC Track
        edm4hep::MutableTrack dcTrack=dcTrackCol->create();

        //Create TrackState
        edm4hep::TrackState trackStateFirstDCHit;
        float charge=trackStateMc.omega/fabs(trackStateMc.omega);
        if(m_useFirstHitForDC&&getTrackStateFirstHit(m_DCSimTrackerHitCol,
                    charge,trackStateFirstDCHit)){
            dcTrack.addToTrackStates(trackStateFirstDCHit);
            dcTrack.addToTrackStates(trackStateMc);
        }else{
            dcTrack.addToTrackStates(trackStateMc);
            dcTrack.addToTrackStates(trackStateFirstDCHit);
        }
        //sdtTk.addToTrackStates(trackStateFirstDCHit);

        ///Add other track properties
        dcTrack.setNdf(dcTrack.trackerHits_size()-5);

        ///Add DC hits to tracks after track state set
        if(m_useIdealHit){
            nDCHitDCTk=addIdealHitsToTk(m_DCDigiCol,truthTrackerHitCol,dcTrack,
                    "DC digi",nDCHitDCTk);
        }else{
            nDCHitDCTk=addHitsToTk(m_DCDigiCol,dcTrack,"DC digi",nDCHitDCTk);
        }
//        if(m_useSi) nDCHitSDTTk=addHitsToTk(m_DCDigiCol,sdtTk,"DC digi",nDCHitSDTTk);
       // if(m_useSi)
       // {
        if(m_smearHits){
            m_nSmearDCTrackHit = smearDCTkhit(m_DCDigiCol,
                    m_SmeartruthTrackerHitCol,m_DCSimTrackerHitCol,
                    w_SimSmearHCol,m_DCHitAssociationCol,
                    w_SmearAssociationCol,m_resX,m_resY,m_resZ);
        }

        if(!m_useNoiseHits)
        {
            //nDCHitSDTTk=addHitsToTk(m_DCDigiCol,sdtTk,"DC digi",nDCHitSDTTk);
            nDCHitSDTTk=addHitsToTk(m_SmeartruthTrackerHitCol,sdtTk,"DC digi",nDCHitSDTTk);
        } else {
            m_nNoiseOnSdtTk = makeNoiseHit(SimVec,Vec,AssoVec,
                    m_SmeartruthTrackerHitCol.get(),w_SmearAssociationCol.get());
            if(debugNoiseHitsCol(Vec))std::cout << " Make noise hits successfully!" << std::endl;
            nDCHitSDTTk=addHitsToTk(w_NoiseHitCol,sdtTk,
                    "DC digi",nDCHitSDTTk);
        } 
        //}

        //track.setType();//TODO
        //track.setChi2(gauss(digiDCHitsCol->size-5(),1));//FIXME
        //track.setDEdx();//TODO

        debug()<<"dcTrack nHit "<<dcTrack.trackerHits_size()<<dcTrack<<endmsg;
    }

    ///Set other track parameters
    //sdtTk.setNdf(sdtTk.trackerHits_size()-5);
    //double radiusOfInnermostHit=1e9;
    //edm4hep::Vector3d digiPos=digiDC.getPosition();
    //double r=sqrt(digiPos.x*digiPos.x+digiPos.y*digiPos.y);
    //if(r<radiusOfInnermostHit) radiusOfInnermostHit=r;

    debug()<<"sdtTk nHit "<<sdtTk.trackerHits_size()<<sdtTk<<endmsg;
    debug()<<"nVXDHit "<<nVXDHit<<" nSITHit "<<nSITHit<<" nSETHit "<<nSETHit
        <<" nFTDHit "<<nFTDHit<<" nDCHitSDTTk "<<nDCHitSDTTk<<endmsg;

    if(m_tuple){
        m_nHitOnSdtTkVXD=nVXDHit;
        m_nHitOnSdtTkSIT=nSITHit;
        m_nHitOnSdtTkSET=nSETHit;
        m_nHitOnSdtTkFTD=nFTDHit;
        m_nHitOnSdtTkDC=nDCHitSDTTk;
        //m_nHitOnSdtTk=sdtTk.trackerHits_size();
        debugEvent();
        StatusCode sc=m_tuple->write();
    }

    return StatusCode::SUCCESS;
}

StatusCode TruthTrackerAlg::finalize()
{
    return GaudiAlgorithm::finalize();
}

void TruthTrackerAlg::getTrackStateFromMcParticle(
        const edm4hep::MCParticleCollection* mcParticleCol,
        edm4hep::TrackState& trackState)
{
    ///Convert MCParticle to DC Track and ReconstructedParticle
    debug()<<"MCParticleCol size="<<mcParticleCol->size()<<endmsg;
    for(auto mcParticle : *mcParticleCol){
        /// skip mcParticleVertex do not have enough associated hits TODO
        ///Vertex
        const edm4hep::Vector3d mcParticleVertex=mcParticle.getVertex();//mm
        edm4hep::Vector3f mcParticleVertexSmeared;//mm
        mcParticleVertexSmeared.x=
            CLHEP::RandGauss::shoot(mcParticleVertex.x,m_resVertexX);
        mcParticleVertexSmeared.y=
            CLHEP::RandGauss::shoot(mcParticleVertex.y,m_resVertexY);
        mcParticleVertexSmeared.z=
            CLHEP::RandGauss::shoot(mcParticleVertex.z,m_resVertexZ);
        ///Momentum
        const edm4hep::Vector3f mcParticleMom=mcParticle.getMomentum();//GeV
        double mcParticlePt=sqrt(mcParticleMom.x*mcParticleMom.x+
                mcParticleMom.y*mcParticleMom.y);
        //double mcParticlePtSmeared=
        //    CLHEP::RandGauss::shoot(mcParticlePt,m_resPT);
        double mcParticleMomPhi=atan2(mcParticleMom.y,mcParticleMom.x);
        double mcParticleMomPhiSmeared=
            CLHEP::RandGauss::shoot(mcParticleMomPhi,m_resMomPhi);
        edm4hep::Vector3f mcParticleMomSmeared;
        mcParticleMomSmeared.x=mcParticlePt*cos(mcParticleMomPhiSmeared);
        mcParticleMomSmeared.y=mcParticlePt*sin(mcParticleMomPhiSmeared);
        mcParticleMomSmeared.z=CLHEP::RandGauss::shoot(mcParticleMom.z,m_resPz);

        ///Converted to Helix
        double B[3]={1e9,1e9,1e9};
        m_dd4hepField.magneticField({0.,0.,0.},B);
        HelixClass helix;
        //float pos[3]={mcParticleVertexSmeared.x,
        //    mcParticleVertexSmeared.y,mcParticleVertexSmeared.z};
        //float mom[3]={mcParticleMomSmeared.x,mcParticleMomSmeared.y,
        //    mcParticleMomSmeared.z};
        ////FIXME DEBUG
        double pos[3]={mcParticleVertex.x,mcParticleVertex.y,mcParticleVertex.z};//mm
        double mom[3]={mcParticleMom.x,mcParticleMom.y,mcParticleMom.z};//mm
        helix.Initialize_VP(pos,mom,mcParticle.getCharge(),B[2]/dd4hep::tesla);
        if(m_tuple) {
            for(int ii=0;ii<3;ii++) {
                m_mcMom[ii]=mom[ii];//GeV
                m_mcPos[ii]=pos[ii];//mm
            }
        }

        ///new Track
        trackState.D0=helix.getD0();
        trackState.phi=helix.getPhi0();
        trackState.omega=helix.getOmega();
        trackState.Z0=helix.getZ0();
        trackState.tanLambda=helix.getTanLambda();
        trackState.referencePoint=helix.getReferencePoint();
        std::array<float,15> covMatrix;
        for(int i=0;i<15;i++){covMatrix[i]=1.;}//FIXME
        trackState.covMatrix=covMatrix;

        getCircleFromPosMom(pos,mom,B[2]/dd4hep::tesla,mcParticle.getCharge(),m_helixRadius,m_helixXC,m_helixYC);

        debug()<<"dd4hep::mm "<<dd4hep::mm<<" dd4hep::cm "<<dd4hep::cm<<endmsg;
        debug()<<"mcParticle "<<mcParticle
            <<" helix radius "<<helix.getRadius()<<" "<<helix.getXC()<<" "
            <<helix.getYC()<<" mm "
            <<" myhelix radius "<<m_helixRadius<<" "<<m_helixXC<<" "
            <<m_helixYC<<" mm "
            <<" momMC "<<mom[0]<<" "<<mom[1]<<" "<<mom[2]<<"GeV"
            <<" posMC "<<pos[0]<<" "<<pos[1]<<" "<<pos[2]<<"mm"
            <<" momPhi "<<mcParticleMomPhi
            <<" mcParticleVertex("<<mcParticleVertex<<")mm "
            <<" mcParticleVertexSmeared("<<mcParticleVertexSmeared<<")mm "
            <<" mcParticleMom("<<mcParticleMom<<")GeV "
            <<" mcParticleMomSmeared("<<mcParticleMomSmeared<<")GeV "
            <<" Bxyz "<<B[0]/dd4hep::tesla<<" "<<B[1]/dd4hep::tesla
            <<" "<<B[2]/dd4hep::tesla<<" tesla"<<endmsg;
    }//end loop over MCParticleCol
}//end of getTrackStateFromMcParticle

bool TruthTrackerAlg::getTrackStateFirstHit(
        DataHandle<edm4hep::SimTrackerHitCollection>& dcSimTrackerHitCol,
        float charge,edm4hep::TrackState& trackState)
{

    const edm4hep::SimTrackerHitCollection* col=nullptr;
    col=dcSimTrackerHitCol.get();
    debug()<<"TruthTrackerAlg::getTrackStateFirstHit"<<endmsg;
    debug()<<"simTrackerHitCol size "<<col->size()<<endmsg;
    float minHitTime=1e9;
    if(nullptr!=col||0==col->size()){
        edm4hep::SimTrackerHit firstHit;
        for(auto dcSimTrackerHit:*col){
            const edm4hep::Vector3f mom=dcSimTrackerHit.getMomentum();
            if(abs(sqrt(mom[0]*mom[0]+mom[1]*mom[1]))<m_momentumCut)continue;//yzhang TEMP skip hits with momentum <0.5GeV/c
            if(dcSimTrackerHit.getTime()<minHitTime) {
                minHitTime=dcSimTrackerHit.getTime();
                firstHit=dcSimTrackerHit;
            }
            //debug()<<"simTrackerHit time "<<dcSimTrackerHit.getTime()
            //    <<" pos "<<dcSimTrackerHit.getPosition()
            //    <<" mom "<<dcSimTrackerHit.getMomentum()<<endmsg;
        }
        const edm4hep::Vector3d pos=firstHit.getPosition();
        const edm4hep::Vector3f mom=firstHit.getMomentum();
        debug()<<"first Hit pos "<<pos<<" mom "<<mom<<" time "<<minHitTime<<endmsg;
        float pos_t[3]={(float)pos[0],(float)pos[1],(float)pos[2]};
        float mom_t[3]={(float)mom[0],(float)mom[1],(float)mom[2]};
        ///Converted to Helix
        double B[3]={1e9,1e9,1e9};
        m_dd4hepField.magneticField({0.,0.,0.},B);
        HelixClass helix;
        helix.Initialize_VP(pos_t,mom_t,charge,B[2]/dd4hep::tesla);
        m_helixRadiusFirst=helix.getRadius();
        m_helixXCFirst=helix.getXC();
        m_helixYCFirst=helix.getYC();

        ///new Track
        trackState.D0=helix.getD0();
        trackState.phi=helix.getPhi0();
        trackState.omega=helix.getOmega();
        trackState.Z0=helix.getZ0();
        trackState.tanLambda=helix.getTanLambda();
        trackState.referencePoint=helix.getReferencePoint();
        std::array<float,15> covMatrix;
        for(int i=0;i<15;i++){covMatrix[i]=100.;}//FIXME
        trackState.covMatrix=covMatrix;
        debug()<<"first hit trackState "<<trackState<<endmsg;
        return true;
    }
    return false;
}//end of getTrackStateFirstHit

void TruthTrackerAlg::debugEvent()
{
    if(m_useSi){
        ///Retrieve silicon Track
        const edm4hep::TrackCollection* siTrackCol=nullptr;
        siTrackCol=m_siSubsetTrackCol.get();
        if(nullptr!=siTrackCol){
            for(auto siTk:*siTrackCol){
                debug()<<"siTk: "<<siTk<<endmsg;
                edm4hep::TrackState trackStat=siTk.getTrackStates(0);//FIXME?
                double B[3]={1e9,1e9,1e9};
                m_dd4hepField.magneticField({0.,0.,0.},B);
                HelixClass helix;
                helix.Initialize_Canonical(trackStat.phi, trackStat.D0, trackStat.Z0,
                        trackStat.omega, trackStat.tanLambda, B[2]/dd4hep::tesla);

                if(m_tuple){
                    m_siMom[0]=helix.getMomentum()[0];
                    m_siMom[1]=helix.getMomentum()[1];
                    m_siMom[2]=helix.getMomentum()[2];
                    m_siPos[0]=helix.getReferencePoint()[0];
                    m_siPos[1]=helix.getReferencePoint()[1];
                    m_siPos[2]=helix.getReferencePoint()[2];
                    m_nHitOnSiTkVXD=nHotsOnTrack(siTk,lcio::ILDDetID::VXD);
                    m_nHitOnSiTkSIT=nHotsOnTrack(siTk,lcio::ILDDetID::SIT);
                    m_nHitOnSiTkSET=nHotsOnTrack(siTk,lcio::ILDDetID::SET);
                    m_nHitOnSiTkFTD=nHotsOnTrack(siTk,lcio::ILDDetID::FTD);
                }
            }//end of loop over siTk
        }
        if(m_tuple){
            //SimTrackerHits
            m_nSimTrackerHitVXD=simTrackerHitColSize(m_VXDCollection);
            m_nSimTrackerHitSIT=simTrackerHitColSize(m_SITCollection);
            m_nSimTrackerHitSET=simTrackerHitColSize(m_SETCollection);
            m_nSimTrackerHitFTD=simTrackerHitColSize(m_FTDCollection);

            //TrackerHits
            m_nTrackerHitVXD=trackerHitColSize(m_VXDTrackerHits);
            m_nTrackerHitSIT=trackerHitColSize(m_SITTrackerHits);
            m_nTrackerHitSET=trackerHitColSize(m_SETTrackerHits);
            m_nTrackerHitFTD=trackerHitColSize(m_FTDTrackerHits);
            m_nTrackerHitDC=trackerHitColSize(m_DCDigiCol);

            //SpacePoints
            if(m_useSiSpacePoint){
                m_nSpacePointSIT=trackerHitColSize(m_SITSpacePointCol);
            }
            //m_nSpacePointSET=trackerHitColSize(m_SETSpacePointCol);
            m_nSpacePointFTD=trackerHitColSize(m_FTDSpacePointCol);
        }
    }
}

int TruthTrackerAlg::addIdealHitsToTk(DataHandle<edm4hep::TrackerHitCollection>&
        colHandle, edm4hep::TrackerHitCollection*& truthTrackerHitCol,
        edm4hep::MutableTrack& track, const char* msg,int nHitAdded)
{
    if(nHitAdded>0) return nHitAdded;
    int nHit=0;
    const edm4hep::TrackerHitCollection* col=colHandle.get();
    debug()<<"add "<<msg<<" "<<col->size()<<" trackerHit"<<endmsg;
    debug()<<track<<endmsg;
    for(auto hit:*col){
        //get end point of this wire
        TVector3 endPointStart(0,0,0);
        TVector3 endPointEnd(0,0,0);
        m_gridDriftChamber->cellposition(hit.getCellID(),endPointStart,
                endPointEnd);//cm

        //calc. doca of helix to wire
        TVector3 wire(endPointStart.X()/dd4hep::mm,endPointStart.Y()/dd4hep::mm,0);//to mm
        TVector3 center(m_helixXC,m_helixYC,0);//mm
        double docaIdeal=(center-wire).Mag()-m_helixRadius;//mm
        TVector3 centerFirst(m_helixXCFirst,m_helixYCFirst,0);//mm
        double docaIdealFirst=(centerFirst-wire).Mag()-m_helixRadiusFirst;//mm

        //add modified hit
        auto tmpHit = truthTrackerHitCol->create();
        tmpHit=hit.clone();
        tmpHit.setTime(fabs(docaIdeal)*1e3/40.);//40#um/ns, drift time in ns
        track.addToTrackerHits(tmpHit);

        long long int detID=hit.getCellID();
        debug()<<" addIdealHitsToTk "<<m_helixRadius<<" center "<<m_helixXC
            <<" "<<m_helixYC<<" mm wire("<<m_decoder->get(detID,"layer")<<","
            <<m_decoder->get(detID,"cellID")<<") "<<wire.X()<<" "<<wire.Y()
            <<"mm docaIdeal "<<docaIdeal<<" docaIdealFirst "<<docaIdealFirst<<"mm "
            <<"hit.Time orignal "<<hit.getTime()<<" new Time "
            <<fabs(docaIdeal)*1e3/40.<<endmsg;
        ++nHit;
    }
    return nHit;
}

bool TruthTrackerAlg::debugNoiseHitsCol(edm4hep::TrackerHitCollection* Vec)
{
    m_nNoiseDCTrackHit = Vec->size();
std::cout << "  m_nNoiseDCTrackHit = " <<  m_nNoiseDCTrackHit << std::endl;
    int ihit=0;
    for(auto Vechit: *Vec)
    {
        m_NoiseDriftDistance[ihit] = (Vechit.getTime())*m_driftVelocity*1e-3;
        ihit++;
    }

    return true;
}

int TruthTrackerAlg::makeNoiseHit(edm4hep::SimTrackerHitCollection* SimVec,
        edm4hep::TrackerHitCollection* Vec,
        edm4hep::MCRecoTrackerAssociationCollection* AssoVec,
        const edm4hep::TrackerHitCollection* digiDCHitsCol,
        const edm4hep::MCRecoTrackerAssociationCollection* assoHits)
{
    int nNoise = 0;
    //loop all hits
    for(unsigned int i = 0; i<(digiDCHitsCol->size()); i++)
    {
        edm4hep::TrackerHit mcHit = digiDCHitsCol->at(i);
        unsigned long long wcellid = mcHit.getCellID();

        float pocaTime = mcHit.getTime();

        // store trackHit
        auto trkHit = Vec->create();

        trkHit.setCellID(wcellid);
        trkHit.setTime(pocaTime);
        trkHit.setEDep(mcHit.getEDep());
        trkHit.setEdx(mcHit.getEdx());
        trkHit.setPosition(mcHit.getPosition());
        trkHit.setCovMatrix(mcHit.getCovMatrix());
        for(int iAsso=0;iAsso<(int) assoHits->size();iAsso++)
        {

            if(assoHits->at(iAsso).getRec().getCellID()==wcellid )
            {
                auto SimtrkHit = SimVec->create();

                SimtrkHit.setCellID(assoHits->at(iAsso).getSim().getCellID());
                SimtrkHit.setEDep(assoHits->at(iAsso).getSim().getEDep());
                SimtrkHit.setTime(assoHits->at(iAsso).getSim().getTime());
                SimtrkHit.setPathLength(assoHits->at(iAsso).getSim().getPathLength());
                SimtrkHit.setQuality(assoHits->at(iAsso).getSim().getQuality());
                SimtrkHit.setPosition(assoHits->at(iAsso).getSim().getPosition());
                SimtrkHit.setMomentum(assoHits->at(iAsso).getSim().getMomentum());
                SimtrkHit.setMCParticle(assoHits->at(iAsso).getSim().getMCParticle());

                auto asso = AssoVec->create();
                asso.setRec(trkHit);
                asso.setSim(SimtrkHit);
                asso.setWeight(SimtrkHit.getEDep()/trkHit.getEDep());
            }
        }

        double prob1 = fRandom.Uniform(1.);
        if(prob1 > m_fHitPurity) continue;

        float noiseTime = m_pocaTime*fRandom.Uniform(1.);
        if(noiseTime>pocaTime) continue;
        nNoise++;
        trkHit.setTime(noiseTime);

        for(int iAsso=0;iAsso<(int) assoHits->size();iAsso++)
        {

            if(assoHits->at(iAsso).getRec().getCellID()==wcellid )
            {

                AssoVec->at(iAsso).setRec(trkHit);

            }

        }

    }

std::cout << " Truth Noise number = " << nNoise << std::endl;
    return nNoise;
}

int TruthTrackerAlg::smearDCTkhit(DataHandle<edm4hep::TrackerHitCollection>&
        colHandle,DataHandle<edm4hep::TrackerHitCollection>& smearCol,
        DataHandle<edm4hep::SimTrackerHitCollection>& SimDCHitCol,
        DataHandle<edm4hep::SimTrackerHitCollection>& SimSmearDCHitCol,
        DataHandle<edm4hep::MCRecoTrackerAssociationCollection>& AssoDCHitCol,
        DataHandle<edm4hep::MCRecoTrackerAssociationCollection>& AssoSmearDCHitCol,
        double resX, double resY, double resZ)
{
    int nHit=0;
    const edm4hep::TrackerHitCollection* col=colHandle.get();
    auto smearTrackerHitCol = smearCol.createAndPut();

    auto simDCCol = SimDCHitCol.get();
    auto SimVec = SimSmearDCHitCol.createAndPut();

    auto assoHits = AssoDCHitCol.get();
    auto AssoVec = AssoSmearDCHitCol.createAndPut();

    //sort,FIXME
    for(auto hit:*col){

        auto smearHit=smearTrackerHitCol->create();

        float truthD = (hit.getTime())*m_driftVelocity*1e-3;
        m_DriftDistance[nHit] = truthD;

        truthD+=gRandom->Gaus(0,fabs(m_resX));
        m_SmearDriftDistance[nHit] = truthD;
        float smearT = (truthD*1e3)/m_driftVelocity;

        smearHit.setTime(smearT);
        smearHit.setCellID(hit.getCellID());
        smearHit.setType(hit.getCellID());
        smearHit.setQuality(hit.getQuality());
        smearHit.setEDep(hit.getEDep());
        smearHit.setEDepError(hit.getEDepError());
        smearHit.setEdx(hit.getEdx());
        smearHit.setPosition(hit.getPosition());
        smearHit.setCovMatrix(hit.getCovMatrix());
        smearHit.addToRawHits(hit.getObjectID());

        ++nHit;
        for(int iAsso=0;iAsso<(int) assoHits->size();iAsso++)
        {

            if(assoHits->at(iAsso).getRec().getCellID()== smearHit.getCellID())
            {
                auto SimtrkHit = SimVec->create();

                SimtrkHit.setCellID(assoHits->at(iAsso).getSim().getCellID());
                SimtrkHit.setEDep(assoHits->at(iAsso).getSim().getEDep());
                SimtrkHit.setTime(assoHits->at(iAsso).getSim().getTime());
                SimtrkHit.setPathLength(assoHits->at(iAsso).getSim().getPathLength());
                SimtrkHit.setQuality(assoHits->at(iAsso).getSim().getQuality());
                SimtrkHit.setPosition(assoHits->at(iAsso).getSim().getPosition());
                SimtrkHit.setMomentum(assoHits->at(iAsso).getSim().getMomentum());
                SimtrkHit.setMCParticle(assoHits->at(iAsso).getSim().getMCParticle());

                auto asso = AssoVec->create();
                asso.setRec(smearHit);
                asso.setSim(SimtrkHit);
                asso.setWeight(SimtrkHit.getEDep()/smearHit.getEDep());
            }
        }
    }
    std::cout << " SimSmear size = " <<  SimVec->size() << " Asso size = " << AssoVec->size() << std::endl;
    std::cout << " Sim size = " <<  simDCCol->size() << " Asso size = " << assoHits->size() << std::endl;

    return nHit;
}

int TruthTrackerAlg::addHitsToTk(DataHandle<edm4hep::TrackerHitCollection>&
        colHandle, edm4hep::MutableTrack& track, const char* msg,int nHitAdded)
{
    if(nHitAdded>0) return nHitAdded;
    int nHit=0;
    const edm4hep::TrackerHitCollection* col=colHandle.get();
    debug()<<"add "<<msg<<" "<<col->size()<<" trackerHit"<<endmsg;
    //sort,FIXME
    for(auto hit:*col){
        track.addToTrackerHits(hit);
        ++nHit;
    }
    return nHit;
}

int TruthTrackerAlg::addSimHitsToTk(
        DataHandle<edm4hep::SimTrackerHitCollection>& colHandle,
        edm4hep::TrackerHitCollection*& truthTrackerHitCol,
        edm4hep::MutableTrack& track, const char* msg,int nHitAdded)
{
    if(nHitAdded>0) return nHitAdded;
    int nHit=0;
    const edm4hep::SimTrackerHitCollection* col=colHandle.get();
    for(auto simTrackerHit:*col){
        auto trackerHit=truthTrackerHitCol->create();
        if(m_skipSecondaryHit&&simTrackerHit.isProducedBySecondary()) {
            debug()<<"skip secondary simTrackerHit "<<msg<<endmsg;
            continue;
        }
        auto& pos = simTrackerHit.getPosition();
        debug()<<" addSimHitsToTk "<<msg<<" "<<sqrt(pos.x*pos.x+pos.y*pos.y)<<endmsg;
        UTIL::BitField64 encoder(lcio::ILDCellID0::encoder_string) ;
        int detID=encoder[lcio::ILDCellID0::subdet] ;
        double resolution[3];
        if(lcio::ILDDetID::VXD==detID){
            for(int i=0;i<3;i++)resolution[i]=m_resVXD[i];
        }else if(lcio::ILDDetID::SIT==detID){
            for(int i=0;i<3;i++)resolution[i]=m_resSIT[i];
        }else if(lcio::ILDDetID::SET==detID){
            for(int i=0;i<3;i++)resolution[i]=m_resSET[i];
        }else if(lcio::ILDDetID::FTD==detID){
            for(int i=0;i<3;i++)resolution[i]=m_resFTDPixel[i];//FIXME
        }else{
            for(int i=0;i<3;i++)resolution[i]=0.003;
        }
        edm4hep::Vector3d posSmeared;//mm
        posSmeared.x=CLHEP::RandGauss::shoot(pos.x,resolution[0]);
        posSmeared.y=CLHEP::RandGauss::shoot(pos.y,resolution[1]);
        posSmeared.z=CLHEP::RandGauss::shoot(pos.z,resolution[2]);
        trackerHit.setPosition(posSmeared) ;
        encoder.setValue(simTrackerHit.getCellID()) ;
        trackerHit.setCellID(encoder.lowWord());//?FIXME
        std::array<float, 6> cov;
        cov[0]=resolution[0]*resolution[0];
        cov[1]=0.;
        cov[2]=resolution[1]*resolution[1];
        cov[3]=0.;
        cov[4]=0.;
        cov[5]=resolution[2]*resolution[2];
        trackerHit.setCovMatrix(cov);
        debug()<<"add simTrackerHit "<<msg<<" trackerHit "<<trackerHit<<endmsg;
        ///Add hit to track
        track.addToTrackerHits(trackerHit);
        trackerHit.setEDep(simTrackerHit.getEDep());
        trackerHit.addToRawHits(simTrackerHit.getObjectID());
        trackerHit.setType(-8);//FIXME?
        ++nHit;
    }
    debug()<<"add simTrackerHit "<<msg<<" "<<nHit<<endmsg;
    return nHit;
}

int TruthTrackerAlg::addHotsToTk(edm4hep::Track& sourceTrack,
        edm4hep::MutableTrack& targetTrack, int hitType,const char* msg,int nHitAdded)
{
    if(nHitAdded>0) return nHitAdded;
    int nHit=0;
    for(unsigned int iHit=0;iHit<sourceTrack.trackerHits_size();iHit++){
        edm4hep::TrackerHit hit=sourceTrack.getTrackerHits(iHit);
        UTIL::BitField64 encoder(lcio::ILDCellID0::encoder_string);
        encoder.setValue(hit.getCellID());
        if(encoder[lcio::ILDCellID0::subdet]==hitType){
            targetTrack.addToTrackerHits(hit);
            debug()<<endmsg<<" add siHit "<<msg<<" "<<iHit<<" "<<hit
                <<" pos "<<hit.getPosition().x<<" "<<hit.getPosition().y<<" "
                <<hit.getPosition().z<<" " <<endmsg;
            ++nHit;
        }
    }
    debug()<<endmsg<<" "<<nHit<<" "<<msg<<" hits add on track"<<endmsg;
    return nHit;
}

int TruthTrackerAlg::nHotsOnTrack(edm4hep::Track& track, int hitType)
{
    int nHit=0;
    for(unsigned int iHit=0;iHit<track.trackerHits_size();iHit++){
        edm4hep::TrackerHit hit=track.getTrackerHits(iHit);
        UTIL::BitField64 encoder(lcio::ILDCellID0::encoder_string);
        encoder.setValue(hit.getCellID());
        if(encoder[lcio::ILDCellID0::subdet]==hitType){
            ++nHit;
        }
    }
    return nHit;
}

int TruthTrackerAlg::trackerHitColSize(DataHandle<edm4hep::TrackerHitCollection>& col)
{
    const edm4hep::TrackerHitCollection* c=col.get();
    if(nullptr!=c) return c->size();
    return 0;
}

int TruthTrackerAlg::simTrackerHitColSize(DataHandle<edm4hep::SimTrackerHitCollection>& col)
{
    const edm4hep::SimTrackerHitCollection* c=col.get();
    if(nullptr!=c) return c->size();
    return 0;
}
//unit length is mm
void TruthTrackerAlg::getCircleFromPosMom(double pos[3],double mom[3],
        double Bz,double q,double& helixRadius,double& helixXC,double& helixYC)
{
    double FCT = 2.99792458E-4;//mm
    double pxy = sqrt(mom[0]*mom[0]+mom[1]*mom[1]);
    helixRadius = pxy / (FCT*Bz);
    double phiMomRefPoint = atan2(mom[1],mom[0]);
    helixXC= pos[0] + helixRadius*cos(phiMomRefPoint-M_PI*0.5*q);
    helixYC= pos[1] + helixRadius*sin(phiMomRefPoint-M_PI*0.5*q);
}
