#include "RecGenfitAlgSDT.h"
#include "GenfitTrack.h"
#include "GenfitFitter.h"
#include "GenfitField.h"

//genfit
#include "EventDisplay.h"

//cepcsw
#include "DetInterface/IGeomSvc.h"
#include "DataHelper/HelixClass.h"
#include "DataHelper/TrackHelper.h"
#include "DetSegmentation/GridDriftChamber.h"
#include "UTIL/ILDConf.h"

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
#include <chrono>
#include "time.h"


DECLARE_COMPONENT( RecGenfitAlgSDT )

/////////////////////////////////////////////////////////////////////
RecGenfitAlgSDT::RecGenfitAlgSDT(const std::string& name,
        ISvcLocator* pSvcLocator):
        GaudiAlgorithm(name, pSvcLocator),m_nPDG(5),m_dd4hep(nullptr),
        m_gridDriftChamber(nullptr),m_decoder(nullptr)
{
    declareProperty("EventHeaderCollection", m_headerCol);
    declareProperty("MCParticleCollection", m_mcParticleCol,
            "Handle of the input MCParticle collection");
    declareProperty("DriftChamberDigiCollection", m_DCDigiCol,
            "Handle of DC digi(TrakerHit) collection");
    declareProperty("DCHitAssociationCollection", m_DCHitAssociationCol,
            "Handle of simTrackerHit and TrackerHit association collection");
    declareProperty("SDTTrackCollection", m_SDTTrackCol,
            "Handle of input silicon track collection");
    declareProperty("SDTRecTrackCollection",m_SDTRecTrackCol,
            "Handle of input silicon rec. track collection");
    declareProperty("DCTrackCollection", m_dcTrackCol,
            "Handle of DC track collection");
    declareProperty("SDTRecParticleCollection", m_SDTRecParticleCol,
            "Handle of silicon+drift chamber rec. particle collection");
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

StatusCode RecGenfitAlgSDT::initialize()
{
    MsgStream log(msgSvc(), name());
    info()<<" RecGenfitAlgSDT initialize()"<<endmsg;
    m_eventNo=0;

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
    m_nRecTrack=0;
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
    NTuplePtr nt(ntupleSvc(), "RecGenfitAlgSDT/recGenfitAlgSDT");
    if(nt){
        m_tuple=nt;
    }else{
      m_tuple=ntupleSvc()->book("RecGenfitAlgSDT/recGenfitAlgSDT",
                CLID_ColumnWiseTuple,"RecGenfitAlgSDT");
        if(m_tuple){
            StatusCode sc;
            sc=m_tuple->addItem("run",m_run);
            sc=m_tuple->addItem("evt",m_evt);
            sc=m_tuple->addItem("tkId",m_tkId);
            sc=m_tuple->addItem("nStdTrack",m_nSdtTrack);

            sc=m_tuple->addItem("nSdtRecTrack",m_nSdtRecTrack);

            sc=m_tuple->addItem("mcIndex",m_mcIndex,0,100);//max. 100 particles
            sc=m_tuple->addItem("seedMomP",m_seedMomP);//for single track debug
            sc=m_tuple->addItem("seedMomPt",m_seedMomPt);
            sc=m_tuple->addItem("seedMomQ",m_seedMomQ);
            sc=m_tuple->addItem("seedMom",3,m_seedMom);
            sc=m_tuple->addItem("seedPos",3,m_seedPos);
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

            sc=m_tuple->addItem("ErrorcovMatrix",15,m_ErrorcovMatrix);
            sc=m_tuple->addItem("D0",m_D0);
            sc=m_tuple->addItem("phi",m_phi);
            sc=m_tuple->addItem("omega",m_omega);
            sc=m_tuple->addItem("Z0",m_Z0);
            sc=m_tuple->addItem("tanLambda",m_tanLambda);

            sc=m_tuple->addItem("mcP_D0",mcP_D0);
            sc=m_tuple->addItem("mcP_phi",mcP_phi);
            sc=m_tuple->addItem("mcP_omega",mcP_omega);
            sc=m_tuple->addItem("mcP_Z0",mcP_Z0);
            sc=m_tuple->addItem("mcP_tanLambda",mcP_tanLambda);

            sc=m_tuple->addItem("pocaPosKal",5,3,m_pocaPosKal);
            sc=m_tuple->addItem("pocaMomKal",5,3,m_pocaMomKal);
            sc=m_tuple->addItem("pocaMomKalP",5,m_pocaMomKalP);
            sc=m_tuple->addItem("pocaMomKalPt",5,m_pocaMomKalPt);
            sc=m_tuple->addItem("chargeKal",5,m_chargeKal);
            sc=m_tuple->addItem("nDofKal",5,m_nDofKal);
            sc=m_tuple->addItem("chi2Kal",5,m_chi2Kal);
            sc=m_tuple->addItem("isFitted",5,m_isFitted);
            sc=m_tuple->addItem("isFitConverged",5,m_isFitConverged);
            sc=m_tuple->addItem("isFitConvergedFully",5,
                    m_isFitConvergedFully);
            sc=m_tuple->addItem("fittedState",5,m_fittedState);
            sc=m_tuple->addItem("nHitFailedKal",5,m_nHitFailedKal);
            sc=m_tuple->addItem("nHitFitted",5,m_nHitFitted);
            sc=m_tuple->addItem("nDCDigi",m_nDCDigi);
            sc=m_tuple->addItem("nHitMc",m_nHitMc);
            sc=m_tuple->addItem("nHitKalInput",m_nHitKalInput,0,30000);
            //10 is greater than # of tracking detectors
            sc=m_tuple->addItem("hitDetID",10,m_nHitDetType);
            sc=m_tuple->addItem("nHitWithFitInfo",5,m_nHitWithFitInfo);
            sc=m_tuple->addItem("nSimDCHit",m_nSimDCHit,0,50000);
            sc=m_tuple->addItem("mdcHitDriftT",m_nSimDCHit,m_mdcHitDriftT);
            sc=m_tuple->addItem("mdcHitDriftDl",m_nSimDCHit,m_mdcHitDriftDl);
            sc=m_tuple->addItem("mdcHitDriftDr",m_nSimDCHit,m_mdcHitDriftDr);
            sc=m_tuple->addItem("mdcHitLr",m_nSimDCHit,m_mdcHitLr);
            sc=m_tuple->addItem("mdcHitLayer",m_nSimDCHit,m_mdcHitLayer);
            sc=m_tuple->addItem("mdcHitWire",m_nSimDCHit,m_mdcHitWire);
            sc=m_tuple->addItem("mdcHitExpDoca",m_nSimDCHit,m_mdcHitExpDoca);
            sc=m_tuple->addItem("mdcHitExpMcDoca",m_nSimDCHit,m_mdcHitExpMcDoca);
            sc=m_tuple->addItem("mdcHitErr",m_nSimDCHit,m_mdcHitErr);
            sc=m_tuple->addItem("exeTime",m_exeTime);
            sc=m_tuple->addItem("mdcHitMcTkId",m_nSimDCHit,m_mdcHitMcTkId);
            sc=m_tuple->addItem("mdcHitMcLr",m_nSimDCHit,m_mdcHitMcLr);
            sc=m_tuple->addItem("mdcHitMcDrift",m_nSimDCHit,m_mdcHitMcDrift);
            sc=m_tuple->addItem("mdcHitMcX",m_nSimDCHit,m_mdcHitMcX);
            sc=m_tuple->addItem("mdcHitMcY",m_nSimDCHit,m_mdcHitMcY);
            sc=m_tuple->addItem("mdcHitMcZ",m_nSimDCHit,m_mdcHitMcZ);
            sc=m_tuple->addItem("mcPocaX",m_nSimDCHit,m_mdcHitExpMcPocaX);
            sc=m_tuple->addItem("mcPocaY",m_nSimDCHit,m_mdcHitExpMcPocaY);
            sc=m_tuple->addItem("mcPocaZ",m_nSimDCHit,m_mdcHitExpMcPocaZ);
            sc=m_tuple->addItem("mcPocaWireX",m_nSimDCHit,m_mdcHitExpMcPocaWireX);
            sc=m_tuple->addItem("mcPocaWireY",m_nSimDCHit,m_mdcHitExpMcPocaWireY);
            sc=m_tuple->addItem("mcPocaWireZ",m_nSimDCHit,m_mdcHitExpMcPocaWireZ);
            debug()<< "Book tuple RecGenfitAlgSDT/recGenfitAlgSDT" << endmsg;
        }else{
            warning()<< "Tuple RecGenfitAlgSDT/recGenfitAlgSDT not booked" << endmsg;
        }
    }//end of book tuple

    //init genfit event display
    //if(m_showDisplay) m_genfitDisplay = genfit::EventDisplay::getInstance();

    return StatusCode::SUCCESS;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

StatusCode RecGenfitAlgSDT::execute()
{
    info()<<"RecGenfitAlgSDT in execute()"<<endmsg;

    std::cout << "m_debug = " << m_debug << std::endl;

    edm4hep::ReconstructedParticleCollection* sdtRecParticleCol=
        m_SDTRecParticleCol.createAndPut();

    edm4hep::TrackCollection* sdtRecTrackCol=
        m_SDTRecTrackCol.createAndPut();

    StatusCode sc=StatusCode::SUCCESS;

    std::cout<<" RecGenfitAlgSDT execute eventNo  "<<m_eventNo++<<std::endl;
    if(m_debug&&(abs(m_eventNoSelection)<1e8)&&m_eventNo!=m_eventNoSelection){
        return sc;
    }

    std::chrono::time_point<std::chrono::high_resolution_clock> start;
    if(m_tuple) start=std::chrono::high_resolution_clock::now();

    /////retrieve EventHeader
    //auto header = _headerCol.get()->at(0);
    //int evtNo = header.getEventNumber();
    //int runNo = header.getRunNumber();
    //info()<<"run "<<header.getEventNumber()
    //  <<" "<<header.getRunNumber()<<std::endl;

    ///retrieve silicon Track and TrackHits
    const edm4hep::TrackCollection* sdtTrackCol=nullptr;
    if(m_SDTTrackCol.exist())sdtTrackCol=m_SDTTrackCol.get();
    if(nullptr==sdtTrackCol || sdtTrackCol->size()<=0) {
        debug()<<"TrackCollection not found or sdtTrackCol size=0"<<endmsg;
        return StatusCode::SUCCESS;
    }

    auto dcHitAssociationCol=m_DCHitAssociationCol.get();
    double eventStartTime=0;

    const edm4hep::TrackCollection* dcTrackCol=nullptr;
    if(m_dcTrackCol.exist()) dcTrackCol=m_dcTrackCol.get();
    if(nullptr==dcTrackCol) {
        debug()<<"TrackCollection not found"<<endmsg;
        return StatusCode::SUCCESS;
    }
    const edm4hep::MCParticleCollection* mcParticleCol=nullptr;
    mcParticleCol=m_mcParticleCol.get();//FIXME get error when call exist()
    if(nullptr==mcParticleCol){
        debug()<<"MCParticleCollection not found"<<endmsg;
        return StatusCode::SUCCESS;
    }
    ///----------------------------------------------------
    ///Loop over Track and do fitting for each track
    ///----------------------------------------------------
    m_firstTuple=true;
    debug()<<"SDTTrackCol size="<<sdtTrackCol->size()<<endmsg;

    std::cout <<"SDTTrackCol size="<<sdtTrackCol->size()<< std::endl;

    for(auto sdtTrack: *sdtTrackCol){
        ///Loop over 5 particle hypothesis(0-4): e,mu,pi,K,p
        ///-1 for chargedgeantino
        for(unsigned int pidType=0;pidType<m_nPDG;pidType++){
            if((m_debugPid>=0) && (m_debugPid!=pidType)) continue;
            debug()<<"processing pidType "<<pidType<<endmsg;
            ///-----------------------------------
            ///Create a GenFit track
            ///-----------------------------------
            GenfitTrack* genfitTrack=new GenfitTrack(m_genfitField,
                    m_gridDriftChamber,m_geomSvc);
            genfitTrack->setDebug(m_debug.value());
            if(m_useTruthTrack){
                //single track only FIXME
                if(!genfitTrack->createGenfitTrackFromMCParticle(pidType,
                            *(mcParticleCol->begin()), eventStartTime)){
                    debug()<<"createGenfitTrackFromMCParticle failed!"<<endmsg;
                    return StatusCode::SUCCESS;
                }
            }else{
                if(!genfitTrack->createGenfitTrackFromEDM4HepTrack(pidType,
                            sdtTrack, eventStartTime,m_isUseCovTrack)){
                    debug()<<"createGenfitTrackFromEDM4HepTrack from SDT track\
                        failed!"<<endmsg;
                    return StatusCode::SUCCESS;
                }
            }
            int nHitAdded=genfitTrack->addHitsOnEdm4HepTrack(sdtTrack,
                    dcHitAssociationCol,m_sigmaHit.value(),
                    m_smearHit.value(),m_fitSiliconOnly.value()
                    ,m_isUseFixedSiHitError.value());
            if(0==nHitAdded){
                debug()<<"No simTrackerHit on track added"<<endmsg;
                return StatusCode::SUCCESS;
            }
            if(m_debug) genfitTrack->printSeed();

            ///-----------------------------------
            ///call genfit fitting procedure
            ///-----------------------------------
            m_genfitFitter->setDebug(m_debug);
            m_genfitFitter->processTrack(genfitTrack,m_resortHits.value());

            ///-----------------------------------
            ///Store track
            ///-----------------------------------
            auto dcRecParticle=sdtRecParticleCol->create();
            auto dcRecTrack=sdtRecTrackCol->create();
            if(!genfitTrack->storeTrack(dcRecParticle,dcRecTrack,pidType,m_ndfCut,
                        m_chi2Cut)){
                debug()<<"Fitting failed!"<<std::endl;
            }else{
                ++m_fitSuccess[pidType];
            }
            if(m_debug) genfitTrack->printSeed();

            if(m_tuple) debugTrack(pidType,genfitTrack);
            if(m_showDisplay) {
                //m_genfitDisplay->addEvent(genfitTrack->getTrack());
                //m_genfitDisplay->open();
            }else{
                delete genfitTrack;
            }
        }//end loop over particle type
    }//end loop over a track
    m_nRecTrack++;

    if(m_tuple) {
        auto finish = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = finish - start;
        info() << "Elapsed time: " << elapsed.count() << " s"<<endmsg;
        m_exeTime=elapsed.count();
        debugEvent(sdtTrackCol,sdtRecTrackCol,eventStartTime);
    }



    //if(m_genfitDisplay) while(1){
    //    std::cout<<"Press any key to finish..."<<std::endl;
    //    //system ("pause");
    //}

    if(m_tuple) sc=m_tuple->write();

    return StatusCode::SUCCESS;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

StatusCode RecGenfitAlgSDT::finalize()
{
    MsgStream log(msgSvc(), name());
    info()<< " RecGenfitAlgSDT in finalize()" << endmsg;

    m_genfitFitter->writeHist();
    delete m_genfitFitter;
    info()<<"RecGenfitAlgSDT nRecTrack="<<m_nRecTrack<<" success e "
        <<m_fitSuccess[0]<<" mu "<<m_fitSuccess[1]<<" pi "<<m_fitSuccess[2]
        <<" K "<<m_fitSuccess[3]<<" p "<<m_fitSuccess[4]<<std::endl;
    if(m_nRecTrack>0){
        std::cout<<"RecGenfitAlgSDT Success rate = "<<std::endl;
        for (int i=0;i<5;i++){
            std::cout<<Form("%d: %d/%d= %2.2f",i,m_fitSuccess[i],m_nRecTrack,
                    ((float) m_fitSuccess[i])/m_nRecTrack)<<std::endl;
        }
    }
    return StatusCode::SUCCESS;
}

void RecGenfitAlgSDT::debugTrack(int pidType,const GenfitTrack* genfitTrack)
{
    /// Get fit status
    const genfit::FitStatus* fitState = genfitTrack->getFitStatus();
    int charge= fitState->getCharge();

    if(m_firstTuple){
        m_nHitKalInput=genfitTrack->getNumPoints();
        debug()<<"m_nHitKalInput "<<m_nHitKalInput<<endmsg;
        //FIXME read from config file
        int detIDs[5]={1,2,3,5,7};//VXD=1,SIT=2,SET=5;FTD=3,
        for(int i=0;i<5;i++){
            m_nHitDetType[detIDs[i]]=genfitTrack->getNumPointsDet(detIDs[i]);
            debug()<<"hot detTypeID "<<detIDs[i]<<" "<<m_nHitDetType[detIDs[i]]<<endmsg;
        }
        m_firstTuple=false;
    }
    m_chargeKal[pidType]= charge;
    m_nHitWithFitInfo[pidType]=genfitTrack->getNumPointsWithFittedInfo();
    m_chi2Kal[pidType]=fitState->getChi2();
    m_nDofKal[pidType]=fitState->getNdf();
    m_isFitted[pidType]=(int)fitState->isFitted();
    m_isFitConverged[pidType]=(int) fitState->isFitConverged();
    m_isFitConvergedFully[pidType]=(int) fitState->isFitConvergedFully();

    ///get fitted state of track
    TMatrixDSym fittedCov;
    TLorentzVector fittedPos;
    TVector3 fittedMom;
    int fittedState=genfitTrack->getFittedState(fittedPos,fittedMom,fittedCov);
    m_fittedState[pidType]=fittedState;
    HelixClass helix;//mm and GeV
    double pos[3]={(fittedPos.X()/dd4hep::mm),(fittedPos.Y()/dd4hep::mm),
        (fittedPos.Z()/dd4hep::mm)};
    double mom[3]={(fittedMom.X()),(fittedMom.Y()),(fittedMom.Z())};
    helix.Initialize_VP(pos,mom,charge,m_genfitField->getBzTesla(fittedPos.Vect()));
    m_pocaMomKalP[pidType]=fittedMom.Mag();

    if(m_debug>0){
        /// Get fit status
        debug()<<"evt "<<m_evt<<" fit result: get status OK? pidType "
            <<pidType<<" fittedState "<<fittedState<<" isFitted "
            <<m_isFitted[pidType]<<" isConverged "<<m_isFitConverged[pidType]
            <<" isFitConvergedFully "<<m_isFitConvergedFully[pidType]
            <<" ndf "<<m_nDofKal[pidType]
            <<" chi2 "<<m_chi2Kal[pidType]<<endmsg;
        if((0!=fittedState)||(!m_isFitted[pidType])||(m_nDofKal[pidType]>m_ndfCut)){
            debug()<<"evt "<<m_evt<<" fit failed"<<endmsg;
        }else{
            debug()<<"evt "<<m_evt<<" fit result: Pos("<<
                fittedPos.X()<<" "<<
                fittedPos.Y()<<" "<<
                fittedPos.Z()<<") mom("<<
                fittedMom.X()<<" "<<
                fittedMom.Y()<<" "<<
                fittedMom.Z()<<") p_tot "<<
                fittedMom.Mag()<<" pt "<<
                fittedMom.Perp()<<endmsg;
        }
    }
}

void RecGenfitAlgSDT::debugEvent(const edm4hep::TrackCollection* sdtTrackCol,
        const edm4hep::TrackCollection* sdtRecTrackCol,
        double eventStartTime)
{
    int iSdtTrack=0;
    m_nSdtTrack=sdtTrackCol->size();
    for(auto sdtTrack: *sdtTrackCol){
        if(iSdtTrack>0) break;//TODO debug for single track only
        edm4hep::TrackState trackStat=sdtTrack.getTrackStates(0);//FIXME?
        HelixClass helixClass;
        helixClass.Initialize_Canonical(trackStat.phi,trackStat.D0,
                trackStat.Z0,trackStat.omega,trackStat.tanLambda,
                m_genfitField->getBzTesla({0.,0.,0.}));

        TLorentzVector posInit(helixClass.getReferencePoint()[0],
                helixClass.getReferencePoint()[1],
                helixClass.getReferencePoint()[2],eventStartTime);
        m_seedPos[0]=posInit.X();
        m_seedPos[1]=posInit.Y();
        m_seedPos[2]=posInit.Z();
        TVector3 momInit(helixClass.getMomentum()[0],
                helixClass.getMomentum()[1],helixClass.getMomentum()[2]);
        m_seedMomP=momInit.Mag();
        m_seedMomPt=momInit.Perp();
        m_seedMom[0]=momInit.X();
        m_seedMom[1]=momInit.Y();
        m_seedMom[2]=momInit.Z();
        iSdtTrack++;
        TVector3 pos,mom;
        TMatrixDSym cov(6);
        double charge;
        CEPC::getPosMomFromTrackState(trackStat,m_genfitField->getBzTesla({0.,0.,0.}),
                pos,mom,charge,cov);
        m_seedMomQ=charge;
        debug()<<" sdtTrack charge "<<charge<<" seed mom "<<momInit.X()<<" "<<
            momInit.Y()<<" "<<momInit.Z()<<endmsg;
        if(m_debug>0){
            pos.Print();
            mom.Print();
            cov.Print();
        }
    }

    const edm4hep::MCParticleCollection* mcParticleCol = nullptr;
    const edm4hep::SimTrackerHitCollection* simDCHitCol=nullptr;

    m_pidIndex=5;

    mcParticleCol=m_mcParticleCol.get();
    int iMcParticle=0;
    HelixClass helix_mcP;
    for(auto mcParticle : *mcParticleCol){
        edm4hep::Vector3f mcPocaMom = mcParticle.getMomentum();//GeV
        edm4hep::Vector3d mcPocaPos = mcParticle.getVertex();

        double mcPos[3]={(mcPocaPos.x),(mcPocaPos.y),(mcPocaPos.z)};
        double mcMom[3]={(mcPocaMom.x),(mcPocaMom.y),(mcPocaMom.z)};
        for(int i=0;i<3;i++){debug()<<"mcPos "<<mcPos[i]<<endmsg;}
        for(int i=0;i<3;i++){debug()<<"mcMom "<<mcMom[i]<<endmsg;}
        float mcCharge = mcParticle.getCharge();
        helix_mcP.Initialize_VP(mcPos,mcMom,mcCharge,
                m_genfitField->getBzTesla(mcPos));

        mcP_D0 = helix_mcP.getD0();
        mcP_phi = helix_mcP.getPhi0();
        mcP_omega = helix_mcP.getOmega();
        mcP_Z0 = helix_mcP.getZ0();
        mcP_tanLambda = helix_mcP.getTanLambda();

        debug()<< " debugEvent Bz " << m_genfitField->getBzTesla(mcPos)
            << " mc d0= " << mcP_D0
            << " phi0= " << mcP_phi
            << " omega= " << mcP_omega
            << " Z0= " << mcP_Z0
            << " tanLambda= " << mcP_tanLambda << endmsg;

        float px=mcPocaMom.x;
        float py=mcPocaMom.y;
        float pz=mcPocaMom.z;
        debug()<<"mc pxyz   "<<px<<" "<<py<<" "<<pz<<endmsg;
        m_pocaMomMcP[iMcParticle]=sqrt(px*px+py*py+pz*pz);
        m_pocaMomMcPt[iMcParticle]=sqrt(px*px+py*py);
        m_pocaMomMc[iMcParticle][0]=px;
        m_pocaMomMc[iMcParticle][1]=py;
        m_pocaMomMc[iMcParticle][2]=pz;
        iMcParticle++;
    }
    m_mcIndex=iMcParticle;

    int iHit=0;
    simDCHitCol=m_simDCHitCol.get();
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
    m_nSimDCHit=simDCHitCol->size();
    const edm4hep::TrackerHitCollection* dCDigiCol=nullptr;
    dCDigiCol=m_DCDigiCol.get();
    if(nullptr!=dCDigiCol){ m_nDCDigi=dCDigiCol->size(); }

    m_nSdtRecTrack=sdtRecTrackCol->size();
    for(auto sdtTrack: *sdtRecTrackCol){
        for(unsigned int i=0; i<sdtTrack.trackStates_size(); i++) {
            edm4hep::TrackState trackStat=sdtTrack.getTrackStates(i);
            std::array<float,15> errorCov;
            errorCov = trackStat.covMatrix;
            for(int j=0; j<15; j++) {
                m_ErrorcovMatrix[j] = errorCov[j];
                debug()<<"debugEvent2 errorCov "<<j<<" "<<errorCov[j]<<endmsg;
            }
            m_D0 = trackStat.D0;
            m_phi = trackStat.phi;
            m_omega = trackStat.omega;
            m_Z0 = trackStat.Z0;
            m_tanLambda = trackStat.tanLambda;

	    debug() << "m_D0" << m_D0 << endmsg; 
	    debug() << "m_phi" << m_phi << endmsg; 
	    debug() << "m_omega" << m_omega << endmsg; 
	    debug() << "m_Z0" << m_Z0 << endmsg; 
	    debug() << "m_tanLambda" << m_tanLambda << endmsg; 
        }
    }
}

void RecGenfitAlgSDT::debugEvent2(const edm4hep::TrackCollection* sdtRecTrackCol) 
{

    m_nSdtRecTrack=sdtRecTrackCol->size();
    for(auto sdtTrack: *sdtRecTrackCol){
        for(int i=0; i<sdtTrack.trackStates_size(); i++) {

            edm4hep::TrackState trackStat=sdtTrack.getTrackStates(i);

            std::array<float,15> errorCov;

            errorCov = trackStat.covMatrix;

            for(int j=0; j<15; j++) {
                m_ErrorcovMatrix[j] = errorCov[j];
            }

             m_D0 = trackStat.D0;
             m_phi = trackStat.phi;
             m_omega = trackStat.omega;
             m_Z0 = trackStat.Z0;
             m_tanLambda = trackStat.tanLambda;
         }

     }
}
