#include "RecGenfitAlgDC.h"
#include "GenfitTrack.h"
#include "GenfitFitter.h"
#include "GenfitField.h"

//genfit
#include "EventDisplay.h"

//cepcsw
#include "DetInterface/IGeomSvc.h"
#include "DataHelper/HelixClass.h"
#include "DetSegmentation/GridDriftChamber.h"

//externals
#include "edm4hep/EventHeaderCollection.h"
#include "edm4hep/MCParticle.h"
#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/SimTrackerHitCollection.h"
#include "edm4hep/TrackerHitCollection.h"
#include "edm4hep/TrackCollection.h"
#include "edm4hep/MCRecoTrackerAssociationCollection.h"
#include "edm4hep/ReconstructedParticle.h"
#include "edm4hep/ReconstructedParticleCollection.h"
#include "edm4hep/Track.h"
#include "edm4hep/TrackCollection.h"
#include "DD4hep/DD4hepUnits.h"
#include "DD4hep/Detector.h"
#include "UTIL/BitField64.h"
#include "DDSegmentation/Segmentation.h"
#include "TRandom.h"
#include "TLorentzVector.h"

//stl
#include "time.h"


DECLARE_COMPONENT( RecGenfitAlgDC )

/////////////////////////////////////////////////////////////////////
RecGenfitAlgDC::RecGenfitAlgDC(const std::string& name, ISvcLocator* pSvcLocator):
        GaudiAlgorithm(name, pSvcLocator),m_nPDG(5),m_dd4hep(nullptr),
        m_gridDriftChamber(nullptr),m_decoder(nullptr)
{
    //declareProperty("EventHeaderCol", _headerCol);
    declareProperty("MCParticle", m_mcParticleCol,
            "Handle of the input MCParticle collection");
    declareProperty("DriftChamberHitsCollection", m_simDCHitCol,
            "Handle of the input SimHit collection");
    declareProperty("DigiDCHitCollection", m_digiDCHitsCol,
            "Handle of digi DCHit collection");
    declareProperty("DCTrackCollection", m_dcTrackCol,
            "Handle of DC track collection");
    //declareProperty("DCHitAssociationCollection", _dcHitAssociationCol,
    //    "Handle of association collection");
    declareProperty("DCRecParticleCollection", m_dcRecParticleCol,
            "Handle of drift chamber reconstructed particle collection");
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

StatusCode RecGenfitAlgDC::initialize()
{
    MsgStream log(msgSvc(), name());
    info()<<" RecGenfitAlgDC initialize()"<<endmsg;

    ///Get GeomSvc
    m_geomSvc=Gaudi::svcLocator()->service("GeomSvc");
    if (nullptr==m_geomSvc) {
        std::cout<<"Failed to find GeomSvc"<<std::endl;
        return StatusCode::FAILURE;
    }
    ///Get Detector
    m_dd4hep = m_geomSvc->lcdd();
    ///Get Field
    m_dd4hepField=m_geomSvc->lcdd()->field();

    /// New a genfit fitter
    m_genfitFitter=new GenfitFitter(m_fitterType.toString().c_str());
    m_genfitField=new GenfitField(m_dd4hepField);
    m_genfitFitter->setField(m_genfitField);
    m_genfitFitter->setGeoMaterial(m_geomSvc->lcdd());
    m_genfitFitter->setEnergyLossBrems(m_correctBremsstrahlung);
    m_genfitFitter->setNoiseBrems(m_correctBremsstrahlung);
    if(m_debug>10) m_genfitFitter->setDebug(m_debug-10);
    if(m_noMaterialEffects) m_genfitFitter->setNoEffects(true);
    if(-1==m_debugPid) m_genfitFitter->setNoEffects(true);
    if(-1==m_debugPid) m_debugPid=0;//charged geantino with electron pid
    if(m_fitterType=="DAF"||m_fitterType=="DafRef"){
        m_genfitFitter->setMaxIterationsBetas(m_bStart,m_bFinal,m_maxIteration);
    } else {
        m_genfitFitter->setMaxIterations(m_maxIteration);
    }
    //print genfit parameters
    if(m_debug) m_genfitFitter->print();
    if(""!=m_genfitHistRootName) m_genfitFitter->initHist(m_genfitHistRootName);

    //initialize member vairables
    for(int i=0;i<5;i++) m_fitSuccess[i]=0;
    m_nDCTrack=0;
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


    ///book tuple
    NTuplePtr nt(ntupleSvc(), "RecGenfitAlgDC/genfit");
    if(nt){
        m_tuple=nt;
    }else{
        m_tuple=ntupleSvc()->book("RecGenfitAlgDC/genfit",
                CLID_ColumnWiseTuple,"RecGenfitAlgDC");
        if(m_tuple){
            StatusCode sc;
            sc=m_tuple->addItem("run",m_run);
            sc=m_tuple->addItem("evt",m_evt);
            sc=m_tuple->addItem("tkId",m_tkId);
            sc=m_tuple->addItem("mcIndex",m_mcIndex,0,100);//max. 100 particles
            sc=m_tuple->addItem("truthPocaMc",m_mcIndex,m_truthPocaMc,3);
            sc=m_tuple->addItem("pocaPosMc",m_mcIndex,m_pocaPosMc,3);
            sc=m_tuple->addItem("pocaMomMc",m_mcIndex,m_pocaMomMc,3);
            sc=m_tuple->addItem("pocaMomMcP",m_mcIndex,m_pocaMomMcP);
            sc=m_tuple->addItem("pocaMomMcPt",m_mcIndex,m_pocaMomMcPt);
            sc=m_tuple->addItem("pocaPosMdc",3,m_pocaPosMdc);
            sc=m_tuple->addItem("pocaMomMdc",3,m_pocaMomMdc);
            sc=m_tuple->addItem("index",m_pidIndex, 0, 5);
            sc=m_tuple->addItem("firstPosKalP",5,3,m_firstPosKal);
            sc=m_tuple->addItem("firstMomKalP",5,m_firstMomKalP);
            sc=m_tuple->addItem("firstMomKalPt",5,m_firstMomKalPt);
            sc=m_tuple->addItem("pocaPosKal",5,3,m_pocaPosKal);
            sc=m_tuple->addItem("pocaMomKal",5,3,m_pocaMomKal);
            sc=m_tuple->addItem("pocaMomKalP",5,m_pocaMomKalP);
            sc=m_tuple->addItem("pocaMomKalPt",5,m_pocaMomKalPt);
            sc=m_tuple->addItem("nDofKal",m_pidIndex,m_nDofKal);
            sc=m_tuple->addItem("chi2Kal",m_pidIndex,m_chi2Kal);
            sc=m_tuple->addItem("convergeKal",m_pidIndex,m_convergeKal);
            sc=m_tuple->addItem("isFittedKal",m_pidIndex,m_isFittedKal);
            sc=m_tuple->addItem("nHitFailedKal",m_pidIndex,m_nHitFailedKal);
            sc=m_tuple->addItem("nHitFitted",m_pidIndex,m_nHitFitted);
            sc=m_tuple->addItem("nDigi",m_nDigi);
            sc=m_tuple->addItem("nHitMc",m_nHitMc);
            sc=m_tuple->addItem("nHitKalInput",m_nHitKalInput);
            sc=m_tuple->addItem("nHitMdc",m_nHitMdc,0,6796);
            sc=m_tuple->addItem("mdcHitDriftT",m_nHitMdc,m_mdcHitDriftT);
            sc=m_tuple->addItem("mdcHitDriftDl",m_nHitMdc,m_mdcHitDriftDl);
            sc=m_tuple->addItem("mdcHitDriftDr",m_nHitMdc,m_mdcHitDriftDr);
            sc=m_tuple->addItem("mdcHitLr",m_nHitMdc,m_mdcHitLr);
            sc=m_tuple->addItem("mdcHitLayer",m_nHitMdc,m_mdcHitLayer);
            sc=m_tuple->addItem("mdcHitWire",m_nHitMdc,m_mdcHitWire);
            sc=m_tuple->addItem("mdcHitExpDoca",m_nHitMdc,m_mdcHitExpDoca);
            sc=m_tuple->addItem("mdcHitExpMcDoca",m_nHitMdc,m_mdcHitExpMcDoca);
            sc=m_tuple->addItem("mdcHitErr",m_nHitMdc,m_mdcHitErr);
            sc=m_tuple->addItem("time",m_pidIndex,m_time);
            sc=m_tuple->addItem("mdcHitMcTkId",m_nHitMdc,m_mdcHitMcTkId);
            sc=m_tuple->addItem("mdcHitMcLr",m_nHitMdc,m_mdcHitMcLr);
            sc=m_tuple->addItem("mdcHitMcDrift",m_nHitMdc,m_mdcHitMcDrift);
            sc=m_tuple->addItem("mdcHitMcX",m_nHitMdc,m_mdcHitMcX);
            sc=m_tuple->addItem("mdcHitMcY",m_nHitMdc,m_mdcHitMcY);
            sc=m_tuple->addItem("mdcHitMcZ",m_nHitMdc,m_mdcHitMcZ);
            sc=m_tuple->addItem("mcPocaX",m_nHitMdc,m_mdcHitExpMcPocaX);
            sc=m_tuple->addItem("mcPocaY",m_nHitMdc,m_mdcHitExpMcPocaY);
            sc=m_tuple->addItem("mcPocaZ",m_nHitMdc,m_mdcHitExpMcPocaZ);
            sc=m_tuple->addItem("mcPocaWireX",m_nHitMdc,m_mdcHitExpMcPocaWireX);
            sc=m_tuple->addItem("mcPocaWireY",m_nHitMdc,m_mdcHitExpMcPocaWireY);
            sc=m_tuple->addItem("mcPocaWireZ",m_nHitMdc,m_mdcHitExpMcPocaWireZ);
            debug()<< "Book tuple RecGenfitAlgDC/genfit" << endmsg;
        }else{
            error()<< "Cannot book tuple RecGenfitAlgDC/genfit" << endmsg;
        }
    }//end of book tuple

    //init genfit event display
    //if(m_showDisplay) m_genfitDisplay = genfit::EventDisplay::getInstance();


    return StatusCode::SUCCESS;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

StatusCode RecGenfitAlgDC::execute()
{
    m_timer=clock();
    info()<<" RecGenfitAlgDC in execute()"<<endmsg;

    /////retrieve EventHeader
    //auto header = _headerCol.get()->at(0);
    //int evtNo = header.getEventNumber();
    //int runNo = header.getRunNumber();

    //std::cout<<"run "<<header.getEventNumber()
    //  <<" "<<header.getRunNumber()<<std::endl;

    ///retrieve Track and TrackHits
    const edm4hep::TrackCollection* dcTrackCol=nullptr;
    dcTrackCol=m_dcTrackCol.get();
    if(nullptr==dcTrackCol) {
        debug()<<"TrackCollection not found"<<endmsg;
        return StatusCode::SUCCESS;
    }
    const edm4hep::TrackerHitCollection* didiDCHitsCol=nullptr;
    didiDCHitsCol=m_digiDCHitsCol.get();
    if(nullptr==didiDCHitsCol) {
        debug()<<"DigiDCHitsCollection not found"<<endmsg;
        return StatusCode::SUCCESS;
    }
    edm4hep::ReconstructedParticleCollection* dcRecParticleCol=
        m_dcRecParticleCol.createAndPut();

    ///----------------------------------------------------
    ///Loop over Track and do fitting for each track
    ///----------------------------------------------------
    debug()<<"DCTrackCol size="<<dcTrackCol->size()<<endmsg;
    for(auto dcTrack: *dcTrackCol){
        ///Loop over 5 particle hypothesis(0-4): e,mu,pi,K,p
        ///, -1 for chargedgeantino
        for(unsigned int pidType=0;pidType<m_nPDG;pidType++){
            //if((m_debugPid>=0) && (m_debugPid!=pidType)) continue;
            ///-----------------------------------
            ///Create a GenFit track
            ///-----------------------------------
            GenfitTrack* genfitTrack=new GenfitTrack(m_genfitField,
                    m_gridDriftChamber);
            genfitTrack->setDebug(m_debug);
            double eventStartTime=0;
            if(!genfitTrack->createGenfitTrackFromEDM4HepTrack(pidType,dcTrack,
                        eventStartTime)){
                debug()<<"createGenfitTrack failed!"<<endmsg;
                return StatusCode::SUCCESS;
            }
            if(m_useTruthHit){
                if(0==genfitTrack->addSpacePointMeasurementOnTrack(dcTrack,m_sigma)){
                    debug()<<"addSpacePointMeasurementOnTrack failed!"<<endmsg;
                    return StatusCode::FAILURE;
                }
            }else{
                if(0==genfitTrack->addWireMeasurementOnTrack(dcTrack,m_sigma)){
                    debug()<<"addWireMeasurementOnTrack failed!"<<endmsg;
                    return StatusCode::SUCCESS;
                }
            }
            if(m_debug) genfitTrack->printSeed();

            ///-----------------------------------
            ///call genfit fitting procedure
            ///-----------------------------------
            m_genfitFitter->processTrack(genfitTrack,m_resortHits);
            m_genfitFitter->setDebug(m_debug);

            ///-----------------------------------
            ///Store track
            ///-----------------------------------
            auto dcRecParticle=dcRecParticleCol->create();
            genfitTrack->storeTrack(dcRecParticle,pidType,m_ndfCut,
                    m_chi2Cut);

            if(m_tuple) debugTrack(pidType,genfitTrack);
            if(m_showDisplay) {
                //m_genfitDisplay->addEvent(genfitTrack->getTrack());
                //m_genfitDisplay->open();
            }else{
                delete genfitTrack;
            }
            ++m_fitSuccess[pidType];
        }//end loop over particle type
    }//end loop over a track

    if(m_tuple) debugEvent();


    //if(m_genfitDisplay) while(1){
    //    std::cout<<"Press any key to finish..."<<std::endl;
    //    //system ("pause");
    //}

    if(m_tuple) m_tuple->write();

    return StatusCode::SUCCESS;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

StatusCode RecGenfitAlgDC::finalize()
{
    MsgStream log(msgSvc(), name());
    info()<< " RecGenfitAlgDC in finalize()" << endmsg;

    m_genfitFitter->writeHist();
    delete m_genfitFitter;
    info()<<"RecGenfitAlgDC nRecTrack="<<m_nDCTrack<<" success e "
        <<m_fitSuccess[0]<<" mu "<<m_fitSuccess[1]<<" pi "<<m_fitSuccess[2]
        <<" K "<<m_fitSuccess[3]<<" p "<<m_fitSuccess[4]<<std::endl;
    if(m_nDCTrack>0){
        std::cout<<"RecGenfitAlgDC Success rate = "<<std::endl;
        for (int i=0;i<5;i++){
            std::cout<<Form("%d %2.2f",i,((float) m_fitSuccess[i])/m_nDCTrack)
                <<std::endl;
        }
    }
    return StatusCode::SUCCESS;
}

void RecGenfitAlgDC::debugTrack(int pidType,const GenfitTrack* genfitTrack)
{
    /// Get fit status
    const genfit::FitStatus* fitState = genfitTrack->getFitStatus();
    double chi2 = fitState->getChi2();
    double ndf = fitState->getNdf();
    double charge = fitState->getCharge();
    int isFitted = fitState->isFitted();
    int isConverged = fitState->isFitConverged();
    int isConvergedFully = fitState->isFitConvergedFully();

    ///get fitted state of track
    TMatrixDSym fittedCov;
    TLorentzVector fittedPos;
    TVector3 fittedMom;
    int fittedState=genfitTrack->getFittedState(fittedPos,fittedMom,fittedCov);
    HelixClass helix;//mm and GeV
    float pos[3]={fittedPos.X()/dd4hep::mm,fittedPos.Y()/dd4hep::mm,fittedPos.Z()/dd4hep::mm};
    float mom[3]={fittedMom.X(),fittedMom.Y(),fittedMom.Z()};
    helix.Initialize_VP(pos,mom,charge,m_genfitField->getBz(fittedPos.Vect()));
    m_pocaMomKalP[pidType]=fittedMom.Mag();

}

void RecGenfitAlgDC::debugEvent()
{
    const edm4hep::MCParticleCollection* mcParticleCol = nullptr;
    const edm4hep::SimTrackerHitCollection* simDCHitCol=nullptr;

    mcParticleCol=m_mcParticleCol.get();
    simDCHitCol=m_simDCHitCol.get();
    m_nHitMdc=simDCHitCol->size();
    int iMcParticle=0;
    int iHit=0;
    for(auto mcParticle : *mcParticleCol){
        for(auto simDCHit: *simDCHitCol){
            edm4hep::Vector3d pos=simDCHit.position();
            TVectorD p(3);
            p[0]=pos.x;//no unit conversion here
            p[1]=pos.y;
            p[2]=pos.z;
            m_mdcHitMcX[iHit]=pos.x;
            m_mdcHitMcY[iHit]=pos.y;
            m_mdcHitMcZ[iHit]=pos.z;
            iHit++;
        }
        edm4hep::Vector3f mcPocaMom = mcParticle.getMomentum();//GeV
        float px=mcPocaMom.x;
        float py=mcPocaMom.y;
        float pz=mcPocaMom.z;
        debug()<<"   "<<px<<" "<<py<<" "<<pz<<endmsg;
        m_pocaMomMcP[iMcParticle]=sqrt(px*px+py*py+pz*pz);
        iMcParticle++;
    }
    m_mcIndex=iHit;

    if(m_tuple) m_tuple->write();
}
