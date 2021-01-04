#include "TruthTrackerAlg.h"
#include "DataHelper/HelixClass.h"
#include "edm4hep/EventHeaderCollection.h"
#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/SimTrackerHitCollection.h"
#include "edm4hep/TrackerHitCollection.h"
#include "edm4hep/TrackCollection.h"
#include "edm4hep/MCRecoTrackerAssociationCollection.h"
#include "edm4hep/MCRecoParticleAssociationCollection.h"
#include "edm4hep/ReconstructedParticleCollection.h"
#include "DDSegmentation/Segmentation.h"
#include "DetSegmentation/GridDriftChamber.h"
#include "DetInterface/IGeomSvc.h"
//#include "edm4hep/SimCalorimeterHitCollection.h"
//#include "edm4hep/CaloHitContributionCollection.h"

#include "DD4hep/Detector.h"
#include "DD4hep/IDDescriptor.h"
#include "DD4hep/Plugins.h"
#include "DD4hep/DD4hepUnits.h"

DECLARE_COMPONENT(TruthTrackerAlg)

TruthTrackerAlg::TruthTrackerAlg(const std::string& name, ISvcLocator* svcLoc)
: GaudiAlgorithm(name,svcLoc),m_dd4hep(nullptr),m_gridDriftChamber(nullptr),
    m_decoder(nullptr)
{
    declareProperty("MCParticle", m_mcParticleCol,
            "Handle of the input MCParticle collection");
    declareProperty("DigiDCHitCollection", m_digiDCHitsCol,
            "Handle of digi DCHit collection");
    declareProperty("DCHitAssociationCollection", m_dcHitAssociationCol,
        "Handle of association collection");
    declareProperty("DCTrackCollection", m_dcTrackCol,
        "Handle of association collection");
    declareProperty("DCRecParticleCollection", m_dcRecParticleCol,
            "Handle of drift chamber reconstructed particle collection");
    declareProperty("DCRecParticleAssociationCollection",
            m_dcRecParticleAssociationCol,
            "Handle of drift chamber reconstructed particle collection");
    declareProperty("debug", m_debug=false);
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
    digiDCHitsCol=m_digiDCHitsCol.get();//FIXME DEBUG
    if(nullptr==digiDCHitsCol){
        debug()<<"TrackerHitCollection not found"<<endmsg;
        //return StatusCode::SUCCESS;
    }
    if((int) digiDCHitsCol->size()>m_maxDCDigiCut) return StatusCode::SUCCESS;

    ///Output Track collection
    edm4hep::TrackCollection* dcTrackCol=m_dcTrackCol.createAndPut();
    ////TODO
    /////Output MCRecoTrackerAssociationCollection collection
    //const edm4hep::MCRecoTrackerAssociationCollection*
    //    mcRecoTrackerAssociationCol=nullptr;
    //if(nullptr==mcRecoTrackerAssociationCol){
    //    log<<MSG::DEBUG<<"MCRecoTrackerAssociationCollection not found"
    //        <<endmsg;
    //    return StatusCode::SUCCESS;
    //}
    //mcRecoTrackerAssociationCol=m_mcRecoParticleAssociation.get();

    ///Convert MCParticle to Track and ReconstructedParticle
    if(m_debug){
        debug()<<"MCParticleCol size="<<mcParticleCol->size()<<endmsg;
    }
    for(auto mcParticle : *mcParticleCol){
        //if(fabs(mcParticle.getCharge()<1e-6) continue;//Skip neutral particles
        edm4hep::Vector3d posV= mcParticle.getVertex();
        edm4hep::Vector3f momV= mcParticle.getMomentum();//GeV
        float pos[3]={(float) posV.x, (float) posV.y, (float) posV.z};
        float mom[3]={momV.x, momV.y, momV.z};
        //FIXME
        //pivotToFirstLayer(mcPocaPos,mcPocaMom,firstLayerPos,firstLayerMom);
        double B[3]={1e9,1e9,1e9};
        m_dd4hepField.magneticField({0.,0.,0.},B);
        if(m_debug)std::cout<<" PDG "<<mcParticle.getPDG()<<" charge "
            << mcParticle.getCharge()
            <<" pos "<<pos[0]<<" "<<pos[1]<<" "<<pos[2]
            <<" mom "<<mom[0]<<" "<<mom[1]<<" "<<mom[2]
            <<" Bxyz "<<B[0]/dd4hep::tesla<<" "<<B[1]/dd4hep::tesla
            <<" "<<B[2]/dd4hep::tesla<<" tesla"<<std::endl;

        HelixClass helix;
        helix.Initialize_VP(pos,mom,(int) mcParticle.getCharge(),B[2]/dd4hep::tesla);
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

        if(m_debug){
            std::cout<<" helix  d0 "<<trackState.D0<<
                " phi "<<trackState.phi<<
                " omega "<<trackState.omega<<
                " z0 "<<trackState.Z0<<
                " tanLambda "<<trackState.tanLambda<<
                " referencePoint "<<trackState.referencePoint<< std::endl;
        }
        //new Track
        edm4hep::Track track = dcTrackCol->create();
        track.addToTrackStates(trackState);
        if(nullptr!=digiDCHitsCol) {
            //track.setType();//TODO
            //track.setChi2(trackState->getChi2());
            track.setNdf(digiDCHitsCol->size()-5);
            //track.setDEdx();//TODO
            //track.setRadiusOfInnermostHit();//TODO
            for(unsigned int i=0; i<digiDCHitsCol->size(); i++ ){
                edm4hep::TrackerHit digiDC=digiDCHitsCol->at(i);
                //if(Sim->MCParti!=current) continue;//TODO
                track.addToTrackerHits(digiDC);
            }
        }
        //TODO
        //new ReconstructedParticle
        //recParticle->setType();
        //dcRecParticle->setEnergy();
        //edm4hep::Vector3f momVec3(helix.getMomentum()[0],
        //        helix.getMomentum()[1],helix.getMomentum()[2]);
        //recParticle->setMomentum(momVec3);
    }

    if(m_debug) std::cout<<"Output DCTrack size="<<dcTrackCol->size()<<std::endl;
    return StatusCode::SUCCESS;
}

StatusCode TruthTrackerAlg::finalize()
{
    return GaudiAlgorithm::finalize();
}
