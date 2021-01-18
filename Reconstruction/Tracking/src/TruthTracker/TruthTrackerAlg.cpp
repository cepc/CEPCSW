#include "TruthTrackerAlg.h"
#include "DataHelper/HelixClass.h"
#include "DDSegmentation/Segmentation.h"
#include "DetSegmentation/GridDriftChamber.h"
#include "DetInterface/IGeomSvc.h"
#include "DD4hep/Detector.h"
#include "DD4hep/IDDescriptor.h"
#include "DD4hep/Plugins.h"
#include "DD4hep/DD4hepUnits.h"

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
    declareProperty("MCParticle", m_mcParticleCol,
            "Handle of the input MCParticle collection");
    declareProperty("DigiDCHitCollection", m_DCDigiCol,
            "Handle of DC digi(TrackerHit) collection");
    declareProperty("DCHitAssociationCollection", m_DCHitAssociationCol,
            "Handle of association collection");
    declareProperty("DCTrackCollection", m_DCTrackCol,
            "Handle of DC track collection");
    declareProperty("SiSubsetTrackCollection", m_siSubsetTrackCol,
            "Handle of silicon subset track collection");
    declareProperty("SDTTrackCollection", m_SDTTrackCol,
            "Handle of SDT track collection");
    declareProperty("DCRecParticleCollection", m_DCRecParticleCol,
            "Handle of drift chamber reconstructed particle collection");
    declareProperty("DCRecParticleAssociationCollection",
            m_DCRecParticleAssociationCol,
            "Handle of drift chamber reconstructed particle collection");
}

StatusCode TruthTrackerAlg::initialize()
{
    ///Get geometry
    m_geomSvc=service<IGeomSvc>("GeomSvc");
    if (!m_geomSvc) {
        error() << "Failed to get GeomSvc." << endmsg;
        return StatusCode::FAILURE;
    }
    ///Get Detector
    m_dd4hep=m_geomSvc->lcdd();
    if (nullptr==m_dd4hep) {
        error() << "Failed to get dd4hep::Detector." << endmsg;
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
        error() << "Failed to get the GridDriftChamber" << endmsg;
        return StatusCode::FAILURE;
    }
    ///Get Decoder
    m_decoder = m_geomSvc->getDecoder(m_readout_name);
    if (nullptr==m_decoder) {
        error() << "Failed to get the decoder" << endmsg;
        return StatusCode::FAILURE;
    }

    return GaudiAlgorithm::initialize();
}

StatusCode TruthTrackerAlg::execute()
{
    ///Retrieve MC particle(s)
    const edm4hep::MCParticleCollection* mcParticleCol=nullptr;
    mcParticleCol=m_mcParticleCol.get();
    if(nullptr==mcParticleCol){
        debug()<<"MCParticleCollection not found"<<endmsg;
        return StatusCode::SUCCESS;
    }
    ///Retrieve DC digi
    const edm4hep::TrackerHitCollection* digiDCHitsCol=nullptr;
    digiDCHitsCol=m_DCDigiCol.get();//FIXME DEBUG
    if(nullptr==digiDCHitsCol){
        debug()<<"TrackerHitCollection not found"<<endmsg;
        //return StatusCode::SUCCESS;//FIXME return when no hits in DC + silicon
    }
    if((int) digiDCHitsCol->size()>m_maxDCDigiCut) return StatusCode::SUCCESS;

    ///Output Track collection
    edm4hep::TrackCollection* dcTrackCol=m_DCTrackCol.createAndPut();
    edm4hep::TrackCollection* sdtTrackCol=m_SDTTrackCol.createAndPut();
    edm4hep::ReconstructedParticleCollection* dcRecParticleCol(nullptr);
    if(m_writeRecParticle){
        dcRecParticleCol=m_DCRecParticleCol.createAndPut();
    }
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

    ///Retrieve silicon Track
    const edm4hep::TrackCollection* siTrackCol=nullptr;
    if(m_siSubsetTrackCol.exist()) siTrackCol=m_siSubsetTrackCol.get();
    if(nullptr!=siTrackCol) {
        ///New SDT track
        for(auto siTrack:*siTrackCol){
            edm4hep::Track sdtTrack=sdtTrackCol->create();
            edm4hep::TrackState sdtTrackState;
            edm4hep::TrackState siTrackStat=siTrack.getTrackStates(0);//FIXME?
            sdtTrackState.location=siTrackStat.location;
            sdtTrackState.D0=siTrackStat.D0;
            sdtTrackState.phi=siTrackStat.phi;
            sdtTrackState.omega=siTrackStat.omega;
            sdtTrackState.Z0=siTrackStat.Z0;
            sdtTrackState.tanLambda=siTrackStat.tanLambda;
            sdtTrackState.referencePoint=siTrackStat.referencePoint;
            for(int k=0;k<15;k++){
                sdtTrackState.covMatrix[k]=siTrackStat.covMatrix[k];
            }
            sdtTrack.addToTrackStates(sdtTrackState);
            sdtTrack.setType(siTrack.getType());
            sdtTrack.setChi2(siTrack.getChi2());
            sdtTrack.setNdf(siTrack.getNdf());
            sdtTrack.setDEdx(siTrack.getDEdx());
            sdtTrack.setDEdxError(siTrack.getDEdxError());
            sdtTrack.setRadiusOfInnermostHit(siTrack.getRadiusOfInnermostHit());
            debug()<<"siTrack trackerHits_size="<<siTrack.trackerHits_size()<<endmsg;
            for(unsigned int iSiTackerHit=0;iSiTackerHit<siTrack.trackerHits_size();
                    iSiTackerHit++){
                sdtTrack.addToTrackerHits(siTrack.getTrackerHits(iSiTackerHit));
            }
            //TODO tracks
            for(auto digiDC:*digiDCHitsCol){
                //if(Sim->MCParti!=current) continue;//TODO
                sdtTrack.addToTrackerHits(digiDC);
            }
        }
    }

    ///Convert MCParticle to Track and ReconstructedParticle
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
        edm4hep::Vector3f mcParticleMom=mcParticle.getMomentum();//GeV
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
        mcParticleMomSmeared.z=
            CLHEP::RandGauss::shoot(mcParticleMom.z,m_resPz);

        ///Converted to Helix
        double B[3]={1e9,1e9,1e9};
        m_dd4hepField.magneticField({0.,0.,0.},B);
        HelixClass helix;
        //float pos[3]={mcParticleVertexSmeared.x,
        //    mcParticleVertexSmeared.y,mcParticleVertexSmeared.z};
        //float mom[3]={mcParticleMomSmeared.x,mcParticleMomSmeared.y,
        //    mcParticleMomSmeared.z};
        ////FIXME DEBUG
        float pos[3]={(float)mcParticleVertex.x,
            (float)mcParticleVertex.y,(float)mcParticleVertex.z};
        float mom[3]={(float)mcParticleMom.x,(float)mcParticleMom.y,
            (float)mcParticleMom.z};
        helix.Initialize_VP(pos,mom,mcParticle.getCharge(),B[2]/dd4hep::tesla);

        ///new Track
        edm4hep::Track dcTrack=dcTrackCol->create();
        edm4hep::TrackState trackState;
        trackState.D0=helix.getD0();
        trackState.phi=helix.getPhi0();
        trackState.omega=helix.getOmega();
        trackState.Z0=helix.getZ0();
        trackState.tanLambda=helix.getTanLambda();
        trackState.referencePoint=helix.getReferencePoint();
        std::array<float,15> covMatrix;
        for(int i=0;i<15;i++){covMatrix[i]=999.;}//FIXME
        trackState.covMatrix=covMatrix;
        dcTrack.addToTrackStates(trackState);
        //dcTrack.setType();//TODO
        //dcTrack.setChi2(gauss(digiDCHitsCol->size-5(),1));//FIXME
        dcTrack.setNdf(digiDCHitsCol->size()-5);
        //dcTrack.setDEdx();//TODO
        //set hits
        double radiusOfInnermostHit=1e9;
        debug()<<"digiDCHitsCol size"<<digiDCHitsCol->size()<<endmsg;
        for(auto digiDC : *digiDCHitsCol){
            //if(Sim->MCParti!=current) continue;//TODO
            edm4hep::Vector3d digiPos=digiDC.getPosition();
            double r=sqrt(digiPos.x*digiPos.x+digiPos.y*digiPos.y);
            if(r<radiusOfInnermostHit) radiusOfInnermostHit=r;
            dcTrack.addToTrackerHits(digiDC);
        }
        dcTrack.setRadiusOfInnermostHit(radiusOfInnermostHit);//TODO

        edm4hep::ReconstructedParticle dcRecParticle;
        if(m_writeRecParticle){
            dcRecParticle=dcRecParticleCol->create();
            ///new ReconstructedParticle
            //dcRecParticle.setType();//TODO
            double mass=mcParticle.getMass();
            double p=sqrt(mcParticleMomSmeared.x*mcParticleMomSmeared.x+
                    mcParticleMomSmeared.y*mcParticleMomSmeared.y+
                    mcParticleMomSmeared.z*mcParticleMomSmeared.z);
            dcRecParticle.setEnergy(sqrt(mass*mass+p*p));
            dcRecParticle.setMomentum(mcParticleMomSmeared);
            dcRecParticle.setReferencePoint(mcParticleVertexSmeared);
            dcRecParticle.setCharge(mcParticle.getCharge());
            dcRecParticle.setMass(mass);
            //dcRecParticle.setGoodnessOfPID();//TODO
            //std::array<float>,10> covMat=?;//TODO
            //dcRecParticle.setCovMatrix(covMat);//TODO
            //dcRecParticle.setStartVertex();//TODO
            edm4hep::ParticleID particleID(0,mcParticle.getPDG(),0,1);//FIXME
            dcRecParticle.setParticleIDUsed(particleID);
            dcRecParticle.addToTracks(dcTrack);
        }//end of write RecParticle

        debug()<<"mcParticle "<<mcParticle
            <<" momPhi "<<mcParticleMomPhi
            <<" mcParticleVertex("<<mcParticleVertex<<")mm "
            <<" mcParticleVertexSmeared("<<mcParticleVertexSmeared<<")mm "
            <<" mcParticleMom("<<mcParticleMom<<")GeV "
            <<" mcParticleMomSmeared("<<mcParticleMomSmeared<<")GeV "
            <<" Bxyz "<<B[0]/dd4hep::tesla<<" "<<B[1]/dd4hep::tesla
            <<" "<<B[2]/dd4hep::tesla<<" tesla"<<endmsg;
        debug()<<"trackState:location,D0,phi,omega,Z0,tanLambda"
            <<",referencePoint,cov"<<std::endl<<trackState<<std::endl;
        debug()<<"dcTrack"<<dcTrack<<endmsg;
        if(m_writeRecParticle) debug()<<"dcRecParticle"<<dcRecParticle<<endmsg;
    }//end loop over MCParticleCol

    debug()<<"Output DCTrack size="<<dcTrackCol->size()<<endmsg;
    return StatusCode::SUCCESS;
}

StatusCode TruthTrackerAlg::finalize()
{
    return GaudiAlgorithm::finalize();
}
