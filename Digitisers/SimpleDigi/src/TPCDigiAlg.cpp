/* -*- Mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */

#include "TPCDigiAlg.h"
#include <iostream>
#include <iomanip>
#include <vector>
#include <map>
#include <cmath>
#include <algorithm>
#include <array>

#include <gsl/gsl_randist.h>

#include "DataHelper/Circle.h"
#include "DataHelper/SimpleHelix.h"
#include "DataHelper/LCCylinder.h"

#include "EventSeeder/IEventSeeder.h"
//#include <IMPL/LCFlagImpl.h>
//#include <IMPL/LCRelationImpl.h>
//#include <IMPL/TrackerHitImpl.h>

//stl exception handler
#include <stdexcept>
#include "CLHEP/Units/SystemOfUnits.h"
#include "voxel.h"

#include "GearSvc/IGearSvc.h"
// STUFF needed for GEAR
#include <gear/GEAR.h>
#include <gear/TPCParameters.h>
#include <gear/PadRowLayout2D.h>
#include <gear/BField.h>
//
#include "UTIL/ILDConf.h"

#define TRKHITNCOVMATRIX 6

//using namespace lcio ;
using namespace std ;

DECLARE_COMPONENT(TPCDigiAlg)

bool compare_phi( Voxel_tpc* a, Voxel_tpc* b)  {
  return ( a->getPhiIndex() < b->getPhiIndex() ) ;
}

bool compare_z( Voxel_tpc* a, Voxel_tpc* b) {
  return ( a->getZIndex() < b->getZIndex() ) ;
}


TPCDigiAlg::TPCDigiAlg(const std::string& name, ISvcLocator* svcLoc)
  : GaudiAlgorithm(name, svcLoc),
    _nEvt(0)
{

  declareProperty("HeaderCol", _headerCol);
  // modify processor description
  //_description = "Produces TPC TrackerHit collection from SimTrackerHit collection, smeared in RPhi and Z. A search is made for adjacent hits on a pad row, if they are closer in z and r-phi than the steering parameters _doubleHitResRPhi (default value 2.0 mm) and _doubleHitResZ (default value 5.0 mm) they are considered to overlap. Clusters of hits smaller than _maxMerge (default value 3) are merged into a single tracker hit, with the position given as the average poision of the hits in phi and in z. Clusters which have _maxMerge hits or more are determined to be identifiable as multiple hits, and are not added to the tracker hit collection. This of course means that good hits caught up in a cluster of background hits will be lossed." ;

  // register steering parameters: name, description, class-variable, default value
  declareProperty("TPCCollection",_padRowHitColHdl,
      "The default pad-row based SimTrackerHit collection");
  //registerInputCollection( LCIO::SIMTRACKERHIT,
  //                        "TPCPadRowHitCollectionName" ,
  //                        "Name of the default pad-row based SimTrackerHit collection"  ,
  //                        _padRowHitColName ,
  //                        std::string("TPCCollection") ) ;

  declareProperty("TPCSpacePointCollection",_spacePointColHdl,
      "Additional space point collection which provides additional guide hits between pad row centers.");
  //registerInputCollection( LCIO::SIMTRACKERHIT,
  //                        "TPCSpacePointCollectionName" ,
  //                        "Name of the additional space point collection which provides additional guide hits between pad row centers."  ,
  //                        _spacePointColName ,
  //                        std::string("TPCSpacePointCollection") ) ;

  declareProperty("TPCLowPtCollection",_lowPtHitsColHdl,
      "The LowPt SimTrackerHit collection Produced by Mokka TPC Driver TPC0X");
  //registerInputCollection( LCIO::SIMTRACKERHIT,
  //                        "TPCLowPtCollectionName" ,
  //                        "Name of the LowPt SimTrackerHit collection Produced by Mokka TPC Driver TPC0X"  ,
  //                        _lowPtHitscolName ,
  //                        std::string("TPCLowPtCollection") ) ;
  //declareProperty("MCParticle", _mcColHdl, "The intput MCParticle collection");

  declareProperty("TPCTrackerHitsCol",_TPCTrackerHitsColHdl,
      "The output TrackerHit collection");
  //registerOutputCollection( LCIO::TRACKERHIT,
  //                         "TPCTrackerHitsCol" ,
  //                         "Name of the Output TrackerHit collection"  ,
  //                         _TPCTrackerHitsCol ,
  //                         std::string("TPCTrackerHits") ) ;

  declareProperty("TPCTrackerHitAssCol", _TPCAssColHdl, "Handle of TrackerHit SimTrackHit relation collection");
  //declareProperty("TPCTrackerHitRelations",_outRelCol,
  //    "The TrackerHit SimTrackHit relation collection");
  //registerOutputCollection(LCIO::LCRELATION,
  //                         "SimTrkHitRelCollection",
  //                         "Name of TrackerHit SimTrackHit relation collection",
  //                         _outRelColName,
  //                         std::string("TPCTrackerHitRelations"));

  declareProperty("UseRawHitsToStoreSimhitPointer",
                  _use_raw_hits_to_store_simhit_pointer=bool(false),
                  "Store the pointer to the SimTrackerHits in RawHits (deprecated) ");
  //registerProcessorParameter("UseRawHitsToStoreSimhitPointer",
  //                           "Store the pointer to the SimTrackerHits in RawHits (deprecated) ",
  //                           _use_raw_hits_to_store_simhit_pointer,
  //                           bool(false));

  declareProperty("PointResolutionPadPhi",_pointResoPadPhi=0.900,
                  "Pad Phi Resolution constant in TPC");
  //registerProcessorParameter( "PointResolutionPadPhi" ,
  //                           "Pad Phi Resolution constant in TPC"  ,
  //                           _pointResoPadPhi ,
  //                           (float)0.900) ;

  declareProperty("RejectCellID0",_rejectCellID0=1,
                  "whether or not to use hits without proper cell ID (pad row)");
  //registerProcessorParameter( "RejectCellID0" ,
  //                           "whether or not to use hits without proper cell ID (pad row)"  ,
  //                           _rejectCellID0 ,
  //                           (int)1) ;

  declareProperty("PointResolutionRPhi",_pointResoRPhi0=0.050,
                  "R-Phi Resolution constant in TPC");
  //registerProcessorParameter( "PointResolutionRPhi" ,
  //                           "R-Phi Resolution constant in TPC"  ,
  //                           _pointResoRPhi0 ,
  //                           (float)0.050) ;

  declareProperty("DiffusionCoeffRPhi",_diffRPhi=0.025,
                  "R-Phi Diffusion Coefficent in TPC");
  //registerProcessorParameter( "DiffusionCoeffRPhi" ,
  //                           "R-Phi Diffusion Coefficent in TPC"  ,
  //                           _diffRPhi ,
  //                           (float)0.025) ;

  declareProperty("N_eff",_nEff=22,
                  "Number of Effective electrons per pad in TPC");
  //registerProcessorParameter( "N_eff" ,
  //                           "Number of Effective electrons per pad in TPC"  ,
  //                           _nEff ,
  //                           (int)22) ;

  declareProperty("PointResolutionZ",_pointResoZ0=0.4,
                  "TPC Z Resolution Coefficent independent of diffusion");
  //registerProcessorParameter( "PointResolutionZ" ,
  //                           "TPC Z Resolution Coefficent independent of diffusion"  ,
  //                           _pointResoZ0 ,
  //                           (float)0.4) ;

  declareProperty("DiffusionCoeffZ",_diffZ=0.08,
                  "Z Diffusion Coefficent in TPC");
  //registerProcessorParameter( "DiffusionCoeffZ" ,
  //                           "Z Diffusion Coefficent in TPC"  ,
  //                           _diffZ ,
  //                           (float)0.08) ;

  declareProperty("HitSortingBinningZ",_binningZ=5.0,
                  "Defines spatial slice in Z");
  //registerProcessorParameter( "HitSortingBinningZ" ,
  //                           "Defines spatial slice in Z"  ,
  //                           _binningZ ,
  //                           (float)5.0) ;

  declareProperty("HitSortingBinningRPhi",_binningRPhi=2.0,
                  "Defines spatial slice in RP");
  //registerProcessorParameter( "HitSortingBinningRPhi" ,
  //                           "Defines spatial slice in RP"  ,
  //                           _binningRPhi ,
  //                           (float)2.0) ;


  declareProperty("DoubleHitResolutionZ",_doubleHitResZ=5.0,
                  "Defines the minimum distance for two seperable hits in Z");
  //registerProcessorParameter( "DoubleHitResolutionZ" ,
  //                           "Defines the minimum distance for two seperable hits in Z"  ,
  //                           _doubleHitResZ ,
  //                           (float)5.0) ;

  declareProperty("DoubleHitResolutionRPhi",_doubleHitResRPhi=2.0,
                  "Defines the minimum distance for two seperable hits in RPhi");
  //registerProcessorParameter( "DoubleHitResolutionRPhi" ,
  //                           "Defines the minimum distance for two seperable hits in RPhi"  ,
  //                           _doubleHitResRPhi ,
  //                           (float)2.0) ;

  declareProperty("MaxClusterSizeForMerge",_maxMerge=3,
                  "Defines the maximum number of adjacent hits which can be merged");
  //registerProcessorParameter( "MaxClusterSizeForMerge" ,
  //                           "Defines the maximum number of adjacent hits which can be merged"  ,
  //                           _maxMerge ,
  //                           (int)3) ;
}


StatusCode TPCDigiAlg::initialize()
{
  debug() << "in TPCDigiAlg initialize()" <<endmsg;

//  // From GNU documentation:
//  // A replacement for the standard terminate_handler which prints
//  // more information about the terminating exception (if any) on stderr. Call ...
//  std::set_terminate (__gnu_cxx::__verbose_terminate_handler);
//
//#ifdef DIGIPLOTS
//  /// Hook an AIDA implementation -----------------------------------------------
//
//  // First create a pointer to the "IAnalysisFactory" of a specific AIDA
//  // implementation. This factory can then be used to produce all other
//  // factories.
//  _AF = AIDA_createAnalysisFactory();
//
//  // Create a ITreeFactory. -----------------------------------------------------
//  // A ITree can be used to store AIDA objects in memory or on disk.
//
//  _TRF = _AF->createTreeFactory();
//
//  /// Create a ITree object which is bound to a file. ---------------------------
//  // You must always create a "ITree" object to create any other factory.
//  /*
//   * Creates a new tree and associates it with a store.
//   * The store is assumed to be read/write.
//   * The store will be created if it does not exist.
//   * @param storeName The name of the store, if empty (""), the tree is
//   *                  created in memory and therefore will not be associated
//   *                  with a file.
//   * @param storeType Implementation specific string, may control store type
//   * @param readOnly If true the store is opened readonly, an exception if it
//   *                 does not exist
//   * @param createNew If false the file must exist, if true the file will be
//   *                  created
//   * @param options Other options, currently are not specified
//   */
//  // ITree * ITreeFactory::create(const std::string & storeName,
//  //                              const std::string & storeType = "",
//  //                              bool readOnly = false,
//  //                              bool createNew = false,
//  //                              const std::string & options = "") ;
//
//  _TREE = _TRF->create("TPCDigi.root",
//                       "root",
//                       false,
//                       true);
//
//  /// Create an IHistogramFactory which is bound to the tree "*_TREE". -----------
//
//  /*
//   * Create an IHistogramFactory.
//   * @param tree The ITree which created histograms will be associated to.
//   * @return     The IHistogramFactory.
//   */
//  // IHistogramFactory * IAnalysisFactory::createHistogramFactory(ITree & tree);
//
//  _HF = _AF->createHistogramFactory(*_TREE);
//
//  _TREE->mkdir("Histograms");
//
//  /*
//   * Create a IHistogram1D.
//   * @param path      The path of the created IHistogram. The path can either
//   *                  be a relative or full path.
//   *                  ("/folder1/folder2/dataName" and
//   *                  "../folder/dataName" are valid paths).
//   *                  All the directories in the path must exist. The
//   *                  characther `/` cannot be used in names; it is only
//   *                  used to delimit directories within paths.
//   * @param title     The title of the IHistogram1D.
//   * @param nBins     The number of bins of the x axis.
//   * @param lowerEdge The lower edge of the x axis.
//   * @param upperEdge The upper edge of the x axis.
//   * @param options   The options for the IHistogram1D. The default is "".
//   *                  "type=efficiency" for an efficiency IHistogram1D.
//   * @return          The newly created IHistogram1D.
//   */
//
//
//
//  _phiDiffHisto = _HF->createHistogram1D("Histograms/phi_diff",
//                                         "Calculated Phi - Track Phi",
//                                         201, -0.05, 0.05);
//
//  _thetaDiffHisto = _HF->createHistogram1D("Histograms/theta_diff",
//                                           "Calculated Theta - Track Theta",
//                                           201, -0.05, 0.05);
//
//  _phiRelHisto = _HF->createHistogram1D("Histograms/padPhi",
//                                        "Phi Relative to the Pad",
//                                        201, 0.0, 6.3);
//
//  _thetaRelHisto = _HF->createHistogram1D("Histograms/padtheta",
//                                          "Theta Relative to the pad",
//                                          201, 0.0, 6.3);
//
//  _rPhiDiffHisto = _HF->createHistogram1D("Histograms/rPhiDiff",
//                                          "rPhi_rec - rPhi_sim",
//                                          201, -1.0, 1.0);
//
//  _zDiffHisto = _HF->createHistogram1D("Histograms/zDiff",
//                                       "Z_rec - Z_sim",
//                                       201, -1.0, 1.0);
//
//  _zPullHisto = _HF->createHistogram1D("Histograms/zPull",
//                                       "(z_rec - z_sim) / Sigma_z",
//                                       201, -10.0, 10.0);
//
//  _phiDistHisto = _HF->createHistogram1D("Histograms/phiDist",
//                                         "phi_rec - Phi_sim",
//                                         201, -1.0, 1.0);
//
//  _rPhiPullHisto = _HF->createHistogram1D("Histograms/rPhiPull",
//                                          "(rPhi_rec - rPhi_sim) / Sigma_rPhi",
//                                          201, -10.0, 10.0);
//
//  _zSigmaVsZHisto = _HF->createHistogram2D("Histograms/zSigmaVsZ",
//                                           "z Sigma vs Z ",
//                                           3000, 0.0, 3000.0,
//                                           201, -0.20, 5.20);
//
//  _zSigmaHisto = _HF->createHistogram1D("Histograms/zSigma",
//                                        "z Sigma ",
//                                        201, -0.20, 5.20);
//
//  _rPhiSigmaHisto = _HF->createHistogram1D("Histograms/rPhiSigma",
//                                           "rPhi Sigma",
//                                           201, -0.20, 0.20);
//
//  _radiusCheckHisto = _HF->createHistogram1D("Histograms/radiusCheck",
//                                             "R_hit - TPC Rmin - ((RowIndex + 0.5 )* padheight)",
//                                             201, -0.20, 0.20);
//
//  _ResidualsRPhiHisto = _HF->createHistogram1D("Histograms/ResidualsRPhi",
//                                               "MC Track Phi - Hit Phi",
//                                               50, -0.001, 0.001);
//
//  _NSimTPCHitsHisto = _HF->createHistogram1D("Histograms/SimTPCHits",
//                                             "Number of SimTPC Hits",
//                                             100, 0.0, 1000000.0);
//
//  _NBackgroundSimTPCHitsHisto = _HF->createHistogram1D("Histograms/NBackgroundSimTPCHits",
//                                                       "Number of Background SimTPC Hits",
//                                                       100, 0.0, 1000000.0);
//
//  _NPhysicsSimTPCHitsHisto = _HF->createHistogram1D("Histograms/NPhysicsSimTPCHits",
//                                                    "Number of PhysicsSimTPC Hits",
//                                                    100, 0.0, 100000.0);
//
//  _NPhysicsAbove02GeVSimTPCHitsHisto = _HF->createHistogram1D("Histograms/NPhysicsAbove02GeVTPCHits",
//                                                              "Number of PhysicsSimTPC Hits above 0.2GeV pt",
//                                                              100, 0.0, 100000.0);
//
//  _NPhysicsAbove1GeVSimTPCHitsHisto = _HF->createHistogram1D("Histograms/NPhysicsAbove1GeVPtTPCHits",
//                                                             "Number of PhysicsSimTPC Hits above 1.0 GeV pt",
//                                                             100, 0.0, 100000.0);
//
//  _NRecTPCHitsHisto = _HF->createHistogram1D("Histograms/NRecTPCHits",
//                                             "Number of Rec TPC Hits",
//                                             50, 0.0, 100000.0);
//
//  _NLostPhysicsTPCHitsHisto = _HF->createHistogram1D("Histograms/NLostPhysicsTPCHits",
//                                                     "Number of PhysicsSimTPC Hits Lost",
//                                                     100, 0.0, 5000.0);
//
//  _NLostPhysicsAbove02GeVPtTPCHitsHisto = _HF->createHistogram1D("Histograms/NLostPhysicsAbove02GeVPtTPCHits",
//                                                                 "Number of PhysicsSimTPC Hits Lost above 0.2 GeV pt",
//                                                                 100, 0.0, 5000.0);
//
//  _NLostPhysicsAbove1GeVPtTPCHitsHisto = _HF->createHistogram1D("Histograms/NLostPhysicsAbove1GeVPtTPCHits",
//                                                                "Number of PhysicsSimTPC Hits Lost above 1.0 GeV pt",
//                                                                100, 0.0, 1000.0);
//
//  _NRevomedHitsHisto = _HF->createHistogram1D("Histograms/NRevomedHits",
//                                              "Number of Removed TPC hits",
//                                              100, 0.0, 1000000.0);
//
//
//  _NKeptPhysicsTPCHitsHistoPercent = _HF->createHistogram1D("Histograms/NKeptPhysicsTPCHitsPercent",
//                                                            "Number of PhysicsSimTPC Hits Kept",
//                                                            303, 0.0, 1.01);
//
//  _NKeptPhysicsAbove02GeVPtTPCHitsHistoPercent = _HF->createHistogram1D("Histograms/NKeptPhysicsAbove02GeVPtTPCHitsPercent",
//                                                                        "Number of PhysicsSimTPC Hits Kept above 0.2 GeV pt",
//                                                                        303, 0.0, 1.01);
//
//  _NKeptPhysicsAbove1GeVPtTPCHitsHistoPercent = _HF->createHistogram1D("Histograms/NKeptPhysicsAbove1GeVPtTPCHitsPercent",
//                                                                       "Number of PhysicsSimTPC Hits Kept above 1.0 GeV pt",
//                                                                       303, 0.0, 1.01);
//
//#endif
//
//
//  printParameters() ;

  // get the GEAR manager
  auto _gear = service<IGearSvc>("GearSvc");
  if ( !_gear ) {
    error() << "Failed to find GearSvc ..." << endmsg;
    return StatusCode::FAILURE;
  }
  _GEAR = _gear->getGearMgr();
  

  // initialize gsl random generator
  _random = gsl_rng_alloc(gsl_rng_ranlxs2);
  _SEEDER = service<IEventSeeder>("EventSeeder");
  _SEEDER->registerAlg(this);

  _cellid_encoder = 0 ;
  _nRun = 0 ;
  _nEvt = 0 ;

  debug()<<" __LINE__"<<endmsg;
  //_cellid_encoder = new BitField64( lcio::ILDCellID0::encoder_string ) ;
  debug()<<" __LINE__"<<endmsg;
  return GaudiAlgorithm::initialize();
}

//TPCDigiAlg::processRunHeader( LCRunHeader* run)
//{
//  _nRun++ ;
//}

StatusCode TPCDigiAlg::execute()
{
  debug() << "in TPCDigiAlg execute()" <<endmsg;

  //auto header = _headerCol.get()->at(0);
  //int evtNo = header.getEventNumber();
  //int runNo = header.getRunNumber();

  unsigned int thisSeed = _SEEDER->getSeed(this, _nEvt, 0);
  gsl_rng_set( _random, thisSeed);

  //debug() << "seed set to " << thisSeed << " for event number "<< evtNo << endmsg;

   int numberOfVoxelsCreated(0);

  _NSimTPCHits = 0;
  _NBackgroundSimTPCHits = 0;
  _NPhysicsSimTPCHits = 0;
  _NPhysicsAbove02GeVSimTPCHits = 0;
  _NPhysicsAbove1GeVSimTPCHits = 0;
  _NRecTPCHits = 0;

  _NLostPhysicsTPCHits = 0;
  _NLostPhysicsAbove02GeVPtTPCHits = 0;
  _NLostPhysicsAbove1GeVPtTPCHits = 0;
  _NRevomedHits = 0;

  static bool firstEvent = true;
  _tpcHitMap.clear();
  _tpcRowHits.clear();

  // fg: make sure we have one message in the log file with run and event number for the DBD production ....
  info() << "  =========  processing event "
         << std::setw(9) << _nEvt/*evtNo*/ << " run "
         << std::setw(9) << 0/*runNo*/
         << "  ========= " << endmsg;


  if(firstEvent==true) {
    if (! _use_raw_hits_to_store_simhit_pointer ) {

      debug() << "The relations to SimTrackerHits are now stored in relation collection TPCTrackerHitRelations\n SimTrackerHits are no longer stored in RawTrackerHits. Enable this deprecated feature by setting UseRawHitsToStoreSimhitPointer to true in steering file." << endmsg;

    }
    else{

        debug() << "SimTrackerHits will be stored in RawTrackerHits. This is a deprecated please use the relations stored in TPCTrackerHitRelations" << endmsg;

    }


  }

  firstEvent = false ;

  const gear::TPCParameters& gearTPC = _GEAR->getTPCParameters() ;
  const gear::PadRowLayout2D& padLayout = gearTPC.getPadLayout() ;
  // this gets the center of the first pad in the pad layout
  const gear::Vector2D padCoord = padLayout.getPadCenter(1) ;

  // this assumes that the pad width in r*phi is the same for all pads
  _padWidth = padLayout.getPadWidth(0)*padCoord[0];
  // set size of row_hits to hold (n_rows) vectors
  _tpcRowHits.resize(padLayout.getNRows());

  //// created the collection which will be written out
  _trkhitVec = _TPCTrackerHitsColHdl.createAndPut();
  _relCol = _TPCAssColHdl.createAndPut();//LCRELATION

  //zhang TODO
  // to store the weights
  //LCFlagImpl lcFlag(0) ;
  //lcFlag.setBit( LCIO::LCREL_WEIGHTED ) ;
  //_relCol->setFlag( lcFlag.getFlag()  ) ;
  //

  //auto mcCol = _mcColHdl.get();
  //auto mcp = mcCol->at(0);
  //debug() << "First MCParticle " << mcp << endmsg;

  // first deal with the pad-row based hits from Mokka
  const edm4hep::SimTrackerHitCollection* STHcol = nullptr;
  try {
    STHcol = _padRowHitColHdl.get();
  }
  catch ( GaudiException &e ) {
    debug() << "Collection " << _padRowHitColHdl.fullKey() << " is unavailable in event " << _nEvt << endmsg;
    return StatusCode::SUCCESS;
  }

  _cellid_encoder = new UTIL::BitField64( lcio::ILDCellID0::encoder_string ) ;
  //int det_id = 0 ;
  //if ( (STHcol != nullptr) && (STHcol->size()>0) ) {
  //  auto SimTHit = STHcol->at( 0 ) ;
  //  _cellid_encoder->setValue(SimTHit.getCellID()) ;
  //  if ( (*_cellid_encoder)[lcio::ILDCellID0::subdet]!= ILDDetID::TPC ){
  //    //fatal() << "unsupported detector ID NOT TPC. det_id = " << det_id << endmsg;
  //    return StatusCode::FAILURE;
  //  }
  //}else{
  //  debug() << "Collection " << _padRowHitCol.fullKey() << " is unavailable in event " << _nEvt << endmsg;
  //  return StatusCode::SUCCESS;
  //}

  if( STHcol != nullptr ){

    int n_sim_hits = STHcol->size()  ;

    //TODO
    //LCFlagImpl colFlag( STHcol->getFlag() ) ;

    _NSimTPCHits = n_sim_hits;

    debug() << "number of Pad-Row based SimHits = " << n_sim_hits << endmsg;
    
    edm4hep::ConstMCParticle nMinus2MCP;
    edm4hep::ConstMCParticle previousMCP;
    edm4hep::ConstSimTrackerHit nMinus2SimHit;
    edm4hep::ConstSimTrackerHit previousSimTHit;

    debug() << "processing nhit=" << n_sim_hits << endmsg;
    // loop over all the pad row based sim hits
    for(unsigned int i = 0; i<n_sim_hits; i++){
      auto SimTHit = STHcol->at(i);

      // this will used for nominaml smearing for very low pt rubish, so set it to zero initially
      double ptSqrdMC = 0;

      float edep;
      double padPhi(0.0);
      double padTheta (0.0);

      auto& pos = SimTHit.getPosition();
      debug() << "processing hit id = " << SimTHit.id() << endmsg;
      debug() << " SimTHit= "
              << " x = "  << pos[0]
              << " y = "  << pos[1]
              << " z = "  << pos[2]
              << endmsg;

      CLHEP::Hep3Vector thisPoint(pos[0],pos[1],pos[2]);
      double padheight = padLayout.getPadHeight(padLayout.getNearestPad(thisPoint.perp(),thisPoint.phi()));

      const double bField = _GEAR->getBField().at( gear::Vector3D( 0., 0., 0.) ).z() ;
      // conversion constant. r = pt / (FCT*bField)
      const double FCT = 2.99792458E-4;
      bool found_mc = false;
      edm4hep::ConstMCParticle mcp;
      try{ // protect crash while MCParticle unavailable
        mcp = SimTHit.getMCParticle() ;
      }
      catch(...){
        debug() << "catch throw MCParticle not available" << endmsg;
      }
      // increase the counters for the different classification of simhits
      //yzhang FIXME mcp!=NULL
      if(mcp.isAvailable()){

        // get the pt of the MCParticle, this will used later to uses nominal smearing for low momentum rubish
        const edm4hep::Vector3f momentumMC= mcp.getMomentum() ;
        ptSqrdMC = momentumMC[0]*momentumMC[0]+momentumMC[1]*momentumMC[1] ;
        
        debug() << " mcp id = " << mcp.id() 
                << " px = "  << momentumMC[0]
                << " py = "  << momentumMC[1]
                << " pz = "  << momentumMC[2]
                << endmsg;
        
        // SJA:FIXME: the fact that it is a physics hit relies on the fact that for overlay
        // the pointer to the mcp is set to NULL. This distinction may not always be true ...
        ++_NPhysicsSimTPCHits ;
        if( ptSqrdMC > (0.2*0.2) ) ++_NPhysicsAbove02GeVSimTPCHits ;
        if( ptSqrdMC > 1.0 )  ++_NPhysicsAbove1GeVSimTPCHits ;
        
#ifdef DIGIPLOTS
//        if(mcp) plotHelixHitResidual(mcp, &thisPoint);
#endif
      } else {
        debug() << "MCParticle not available" << endmsg;
        ++_NBackgroundSimTPCHits;
      }

      // if the hits contain the momentum of the particle use this to calculate the angles relative to the pad
      //yzhang TODO
      //colFlag.bitSet(LCIO::THBIT_MOMENTUM) == true
      //if(colFlag.bitSet(LCIO::THBIT_MOMENTUM))
      
      bool colFlag_bitSet = true;
      if(colFlag_bitSet){
        // FIXME: fucd, since momentum from SimTrackerHit available, change MCParticle's to SimTrackerHit's, which better?
        //const edm4hep::Vector3f mcpMomentum = mcp.getMomentum() ;
        const edm4hep::Vector3f mcpMomentum = SimTHit.getMomentum();

        CLHEP::Hep3Vector mom(mcpMomentum[0],mcpMomentum[1],mcpMomentum[2]);
        // const double pt = mom.perp();
        // const double radius = pt / (FCT*bField);
        
        // const double tanLambda = mom.z()/pt;
        padPhi = fabs(thisPoint.deltaPhi(mom));
        padTheta = mom.theta();
      }
      else { // LCIO::THBIT_MOMENTUM not set

        // as the momentum vector is not available from the hits use triplets of
        // hits to fit a circle and calculate theta and phi relative to the pad

        if ((!mcp.isAvailable()) || (sqrt(ptSqrdMC) / (FCT*bField)) < ( padheight / (0.1 * CLHEP::twopi))) {
          // if the hit has no record of it MCParticle then there is no way to know if this hit has consecutive hits from the same MCParticle
          // so just set nominal values theta=phi=90
          // here make a cut for particles which will suffer more than a 10 percent change in phi over the distance of the pad
          // R > padheight/(0.1*2PI)
          // in both cases set the angles to 90 degrees
          padTheta = CLHEP::twopi/4.0 ;
          padPhi = CLHEP::twopi/4.0 ;
        }
        else{
          edm4hep::ConstSimTrackerHit nextSimTHit;
          edm4hep::ConstSimTrackerHit nPlus2SimHit;
          edm4hep::ConstMCParticle nextMCP;
          edm4hep::ConstMCParticle nPlus2MCP;
          // if there is at least one more hit after this one, set the pointer to the MCParticle for the next hit
          if (i < (n_sim_hits-1) ) {
            nextSimTHit = STHcol->at( i+1 ) ;
            nextMCP     = nextSimTHit.getMCParticle() ;
          }
          else{ // set make sure that the pointers are set back to NULL so that the comparisons later hold
            //nextSimTHit = edm4hep::ConstSimTrackerHit;
            //nextMCP     = edm4hep::ConstMCParticle;
          }
          // if there is at least two more hits after this one, set the pointer to the MCParticle for the next but one hit
          if (i < (n_sim_hits-2) ) {
            nPlus2SimHit = STHcol->at( i+2 );
            nPlus2MCP    = nPlus2SimHit.getMCParticle() ;
          }
          else{ // set make sure that the pointers are set back to NULL so that the comparisons later hold
            //_nPlus2SimHit = edm4hep::ConstSimTrackerHit;
            //_nPlus2MCP    = edm4hep::ConstMCParticle;
          }

          if      ( mcp==previousMCP && mcp==nextMCP )    { // middle hit of 3 from the same MCParticle

            CLHEP::Hep3Vector precedingPoint(previousSimTHit.getPosition()[0],previousSimTHit.getPosition()[1],previousSimTHit.getPosition()[2]) ;
            CLHEP::Hep3Vector followingPoint(nextSimTHit.getPosition()[0],nextSimTHit.getPosition()[1],nextSimTHit.getPosition()[2]) ;

            debug() << "address of _previousSimTHit = " << previousSimTHit
                    << " x = "  << previousSimTHit.getPosition()[0]
                    << " y = "  << previousSimTHit.getPosition()[1]
                    << " z = "  << previousSimTHit.getPosition()[2]
                    << std::endl ;

            debug() << "address of _nextSimTHit = " << nextSimTHit
                    << " x = "  << nextSimTHit.getPosition()[0]
                    << " y = "  << nextSimTHit.getPosition()[1]
                    << " z = "  << nextSimTHit.getPosition()[2]
                    << std::endl ;

            // get phi and theta using functions defined below
            padPhi = getPadPhi( &thisPoint, &precedingPoint, &thisPoint, &followingPoint);
            padTheta = getPadTheta(&precedingPoint, &thisPoint, &followingPoint);

          }
          else if ( mcp==nextMCP     && mcp==nPlus2MCP )  { // first  hit of 3 from the same MCParticle

            CLHEP::Hep3Vector followingPoint(nextSimTHit.getPosition()[0],nextSimTHit.getPosition()[1],nextSimTHit.getPosition()[2]) ;
            CLHEP::Hep3Vector nPlus2Point(nPlus2SimHit.getPosition()[0],nPlus2SimHit.getPosition()[1],nPlus2SimHit.getPosition()[2]) ;

            // get phi and theta using functions defined below
            padPhi = getPadPhi( &thisPoint, &thisPoint, &followingPoint, &nPlus2Point);
            padTheta = getPadTheta(&thisPoint, &followingPoint, &nPlus2Point);

          }
          else if ( mcp==previousMCP && mcp==nMinus2MCP ) { // last   hit of 3 from the same MCParticle
            auto& posMinus = nMinus2SimHit.getPosition();
            auto& posPrevious = previousSimTHit.getPosition();
            
            CLHEP::Hep3Vector nMinus2Point(posMinus[0],posMinus[1],posMinus[2]);
            CLHEP::Hep3Vector precedingPoint(posPrevious[0],posPrevious[1],posPrevious[2]);

            // get phi and theta using functions defined below
            padPhi = getPadPhi( &thisPoint, &nMinus2Point, &precedingPoint, &thisPoint);
            padTheta = getPadTheta(&nMinus2Point, &precedingPoint, &thisPoint);

          }
          else{ // the hit is isolated as either a single hit, or a pair of hits, from a single MCParticle
            padTheta = CLHEP::twopi/4.0 ;
            padPhi = CLHEP::twopi/4.0 ;
          }
        }

//#ifdef DIGIPLOTS
//        if(colFlag.bitSet(LCIO::THBIT_MOMENTUM)) {
//
//          const float * mcpMomentum = SimTHit.getMomentum() ;
//
//          CLHEP::Hep3Vector mom(mcpMomentum[0],mcpMomentum[1],mcpMomentum[2]);
//
//          double trackPhi = mom.phi();
//
//          if(trackPhi<0.0) trackPhi=trackPhi+twopi;
//          if(trackPhi>twopi) trackPhi=trackPhi-twopi;
//          if(trackPhi>twopi/2.0) trackPhi = trackPhi - twopi/2.0 ;
//
//          double localPhi = thisPoint.phi() - padPhi;
//
//          _phiRelHisto->fill(padPhi);
//          _phiDiffHisto->fill((fabs(localPhi - trackPhi))/trackPhi);
//          _thetaRelHisto->fill(padTheta);
//          _thetaDiffHisto->fill( (sin(padTheta) - sin(mom.theta()))/sin(mom.theta()) );
//
//          debug() << "track Phi = " << trackPhi * (360.0 / twopi) << endmsg;
//          debug() << "localPhi = " << localPhi * (360.0 / twopi) << endmsg;
//          debug() << "pad Phi = " << padPhi * (360.0 / twopi) << endmsg;
//          debug() << "pad Phi from track mom = " << ( thisPoint.phi() - trackPhi ) * (360.0 / twopi) << endmsg;
//          debug() << "padTheta = " << padTheta * (360.0 / twopi) << endmsg;
//          debug() << "padTheta from track mom = " << mom.theta() * (360.0 / twopi) << endmsg;
//
//        }
//#endif

      }
      debug() << "padPhi = " << padPhi << " padTheta = " << padTheta << endmsg;
      //      int pad = padLayout.getNearestPad(thisPoint.perp(),thisPoint.phi());
      int layerNumber = SimTHit.getCellID();

      if(_rejectCellID0 && (layerNumber<1)) {
        continue;
      }

      edep = SimTHit.getEDep();

      // Calculate Point Resolutions according to Ron's Formula

      // sigma_{RPhi}^2 = sigma_0^2 + Cd^2/N_{eff} * L_{drift}

      // sigma_0^2 = (50micron)^2 + (900micron*sin(phi))^2
      // Cd^2/N_{eff}} = 25^2/(22/sin(theta)*h/6mm)
      // Cd = 25 ( microns / cm^(1/2) )
      // (this is for B=4T, h is the pad height = pad-row pitch in mm,
      // theta is the polar angle)

      // sigma_{z}^2 = (400microns)^2 + L_{drift}cm * (80micron/sqrt(cm))^2

      double aReso =_pointResoRPhi0*_pointResoRPhi0 + (_pointResoPadPhi*_pointResoPadPhi * sin(padPhi)*sin(padPhi)) ;
      double driftLength = gearTPC.getMaxDriftLength() - (fabs(thisPoint.z()));

      if (driftLength <0) {
        debug() << " TPCDigiAlg : Warning! driftLength < 0 " << driftLength << " --> Check out your GEAR file!!!!" << endmsg;
        debug() << "Setting driftLength to 0.1" << endmsg;
        debug() << "gearTPC.getMaxDriftLength() = " << gearTPC.getMaxDriftLength() << endmsg;
        driftLength = 0.10;
      }

      padheight = padLayout.getPadHeight(padLayout.getNearestPad(thisPoint.perp(),thisPoint.phi()));

      double bReso = ( (_diffRPhi * _diffRPhi) / _nEff ) * sin(padTheta) * ( 6.0 / (padheight) )  * ( 4.0 / bField  ) ;

      double tpcRPhiRes = sqrt( aReso + bReso * (driftLength / 10.0) ); // driftLength in cm

      double tpcZRes  = sqrt(( _pointResoZ0 * _pointResoZ0 )
                             +
                             ( _diffZ * _diffZ ) * (driftLength / 10.0) ); // driftLength in cm

      int padIndex = padLayout.getNearestPad(thisPoint.perp(),thisPoint.phi());

      const gear::DoubleVec & planeExt = padLayout.getPlaneExtent() ;
      double TPCPadPlaneRMin = planeExt[0] ;
      double TPCPadPlaneRMax = planeExt[1] ;

      int iRowHit = padLayout.getRowNumber(padIndex);
      int iPhiHit = padLayout.getPadNumber(padIndex);
      int NBinsZ =  (int) ((2.0 * gearTPC.getMaxDriftLength()) / _binningZ);
      int iZHit = (int) ( (float) NBinsZ * ( gearTPC.getMaxDriftLength() + thisPoint.z() ) / ( 2.0 * gearTPC.getMaxDriftLength() ) ) ;

      if(iZHit<0) iZHit=0;
      if(iZHit>NBinsZ) iZHit=NBinsZ;

      // make sure that the hit lies at the middle of the pad ring
      thisPoint.setPerp(padLayout.getPadCenter(padIndex)[0]);

      if( (thisPoint.perp() < TPCPadPlaneRMin) || (thisPoint.perp() > TPCPadPlaneRMax) ) {
        debug() << "Hit R not in TPC " << endmsg;
        debug() << "R = " << thisPoint.perp() << endmsg;
        debug() << "the tpc InnerRadius = " << TPCPadPlaneRMin << endmsg;
        debug() << "the tpc OuterRadius = " << TPCPadPlaneRMax << endmsg;
        debug() << "Hit Dropped " << endmsg;
        continue;
      }

      if( (fabs(thisPoint.z()) > gearTPC.getMaxDriftLength()) ) {
        debug() << "Hit Z not in TPC " << endmsg;
        debug() << "Z = " << thisPoint.z() << endmsg;
        debug() << "the tpc Max Z = " << gearTPC.getMaxDriftLength() << endmsg;
        debug() << "Hit Dropped " << endmsg;
        continue;
      }

      // create a tpc voxel hit and store it for this row
      Voxel_tpc * atpcVoxel = new Voxel_tpc(iRowHit,iPhiHit,iZHit, thisPoint, edep, tpcRPhiRes, tpcZRes);

      _tpcRowHits.at(iRowHit).push_back(atpcVoxel);
      ++numberOfVoxelsCreated;

      // store the simhit pointer for this tpcvoxel hit in the hit map
      _tpcHitMap[atpcVoxel] = SimTHit;

      // move the pointers on
      nMinus2MCP = previousMCP;
      previousMCP = mcp ;
      nMinus2SimHit = previousSimTHit;
      previousSimTHit = SimTHit;

    }
  }

  // now process the LowPt collection
  try {
    auto STHcolLowPt = _lowPtHitsColHdl.get();

    if(nullptr!=STHcolLowPt){

      int n_sim_hitsLowPt = STHcolLowPt->size()  ;

      _NBackgroundSimTPCHits += n_sim_hitsLowPt;
      _NSimTPCHits += n_sim_hitsLowPt;

      _padWidth = padLayout.getPadWidth(0)*padCoord[0];

      debug() << "number of LowPt hits:" << n_sim_hitsLowPt << endmsg;

      // loop over the LowPt hit collection
      for(auto SimTHit : *STHcolLowPt){

        CLHEP::Hep3Vector thisPoint(SimTHit.getPosition()[0],SimTHit.getPosition()[1],SimTHit.getPosition()[2]);

        const gear::DoubleVec & planeExt = padLayout.getPlaneExtent() ;
        double TPCPadPlaneRMin = planeExt[0] ;
        double TPCPadPlaneRMax = planeExt[1] ;

        int NBinsZ =  (int) ((2.0 * gearTPC.getMaxDriftLength()) / _binningZ);

        if( (thisPoint.perp() < TPCPadPlaneRMin) || (thisPoint.perp() > TPCPadPlaneRMax) ) {
          debug() << "Hit R not in TPC " << endmsg;
          debug() << "R = " << thisPoint.perp() << endmsg;
          debug() << "the tpc InnerRadius = " << TPCPadPlaneRMin << endmsg;
          debug() << "the tpc OuterRadius = " << TPCPadPlaneRMax << endmsg;
          debug() << "Hit Dropped " << endmsg;
          continue;
        }

        if( (fabs(thisPoint.z()) > gearTPC.getMaxDriftLength()) ) {
          debug() << "Hit Z not in TPC " << endmsg;
          debug() << "Z = " << thisPoint.z() << endmsg;
          debug() << "the tpc Max Z = " << gearTPC.getMaxDriftLength() << endmsg;
          debug() << "Hit Dropped " << endmsg;
          continue;
        }

        int padIndex = padLayout.getNearestPad(thisPoint.perp(),thisPoint.phi());
        //double padheight = padLayout.getPadHeight(padIndex);

        int iRowHit = padLayout.getRowNumber(padIndex);
        int iPhiHit = padLayout.getPadNumber(padIndex);
        int iZHit = (int) ( (float) NBinsZ *
            ( gearTPC.getMaxDriftLength() + thisPoint.z() ) / ( 2.0 * gearTPC.getMaxDriftLength() ) ) ;

        // shift the hit in r-phi to the nearest pad-row centre
        thisPoint.setPerp(padLayout.getPadCenter(padIndex)[0]);

        // set the resolutions to the pads to digital like values
        double tpcRPhiRes = _padWidth;
        double tpcZRes = _binningZ;

        // create a tpc voxel hit for this simhit and store it for this tpc pad row
        Voxel_tpc * atpcVoxel = new Voxel_tpc(iRowHit,iPhiHit,iZHit, thisPoint, SimTHit.getEDep(), tpcRPhiRes, tpcZRes);

        _tpcRowHits.at(iRowHit).push_back(atpcVoxel);
        ++numberOfVoxelsCreated;

        // store the simhit pointer for this voxel hit in a map
        _tpcHitMap[atpcVoxel] = SimTHit;

      }
    }
  }catch(GaudiException e){
    warning() << "Catch exception when read LowPtHitsCol" << endmsg;
  }

  int number_of_adjacent_hits(0);

  debug() << "finished looping over simhits, number of voxels = " << numberOfVoxelsCreated << endmsg;

  int numberOfhitsTreated(0);

  vector <Voxel_tpc *> row_hits;

  // loop over the tpc rows containing hits and check for merged hits
  for (unsigned int i = 0; i<_tpcRowHits.size(); ++i){

    row_hits = _tpcRowHits.at(i);
    std::sort(row_hits.begin(), row_hits.end(), compare_phi );

    // double loop over the hits in this row
    for (unsigned int j = 0; j<row_hits.size(); ++j){

      ++numberOfhitsTreated;

      for (unsigned int k = j+1; k<row_hits.size(); ++k){

        if(row_hits[k]->getPhiIndex() > (row_hits[j]->getPhiIndex())+2){ // SJA:FIXME: here we need an OR to catch the wrap around
          break; // only compare hits in adjacent phi bins
        }

        // look to see if the two hit occupy the same pad in phi or if not whether they are within the r-phi double hit resolution
        else if( row_hits[k]->getPhiIndex()==row_hits[j]->getPhiIndex()
                ||
                ( (fabs(row_hits[k]->getHep3Vector().deltaPhi(row_hits[j]->getHep3Vector()))) * row_hits[j]->getR()) < _doubleHitResRPhi ) {

          // if neighboring in phi then compare z
          map <Voxel_tpc*,edm4hep::SimTrackerHit> ::iterator it;

          edm4hep::SimTrackerHit Hit1;
          edm4hep::SimTrackerHit Hit2;

          // search of the simhit pointers in the tpchit map
          it=_tpcHitMap.find(row_hits[j]);
          if(it!= _tpcHitMap.end()) {
            Hit1 = it->second ; // hit found
          }

          it=_tpcHitMap.find(row_hits[k]);
          if(it!= _tpcHitMap.end()) {
            Hit2 = it->second ; // hit found
          }

          double pathlengthZ1(0.0);
          double pathlengthZ2(0.0);

          if( Hit1.isAvailable() && Hit2.isAvailable() ){ // if both sim hits were found

            // check if the track momentum has been stored for the hits
            bool momentum_set = true;

            if( STHcol != NULL ){
              //yzhang TODO
              //LCFlagImpl colFlag( STHcol->getFlag() ) ;
              //momentum_set = momentum_set && colFlag.bitSet(LCIO::THBIT_MOMENTUM) ;
            }

            //if( STHcolLowPt != NULL ){
              //yzhang TODO
              //LCFlagImpl colFlag( STHcolLowPt->getFlag() ) ;
              //momentum_set =  momentum_set && colFlag.bitSet(LCIO::THBIT_MOMENTUM) ;
            //}

            if( momentum_set ){

              const edm4hep::Vector3f Momentum1 = Hit1.getMomentum() ;
              const edm4hep::Vector3f Momentum2 = Hit1.getMomentum() ;

              CLHEP::Hep3Vector mom1(Momentum1[0],Momentum1[1],Momentum1[2]);
              CLHEP::Hep3Vector mom2(Momentum2[0],Momentum2[1],Momentum2[2]);

              pathlengthZ1 = fabs( Hit1.getPathLength() * mom1.cosTheta() );
              pathlengthZ2 = fabs( Hit2.getPathLength() * mom2.cosTheta() );
            }
            else {
              pathlengthZ1 = _doubleHitResZ ; // assume the worst i.e. that the track is moving in z
              pathlengthZ2 = _doubleHitResZ ; // assume the worst i.e. that the track is moving in z
            }

            double dZ = fabs(row_hits[j]->getZ() - row_hits[k]->getZ());

            double spacial_coverage = 0.5*(pathlengthZ1 + pathlengthZ2) + _binningZ;

            if( (dZ - spacial_coverage) < _doubleHitResZ ){

              row_hits[j]->setAdjacent(row_hits[k]);
              row_hits[k]->setAdjacent(row_hits[j]);
              ++number_of_adjacent_hits;

            }
          } else {
            debug() << "Hit1=" << Hit1 << "Hit2=" << Hit2 << endmsg;
          }
        }
      }
    }


    // now all hits have been checked for adjacent hits, go throught and write out the hits or merge

    for (unsigned int j = 0; j<row_hits.size(); ++j){

      Voxel_tpc* seed_hit = row_hits[j];

      if(seed_hit->IsMerged() || seed_hit->IsClusterHit()) {
        continue;
      }

      if(seed_hit->getNumberOfAdjacent()==0){ // no adjacent hits so smear and write to hit collection
        writeVoxelToHit(seed_hit);
      } else if(seed_hit->getNumberOfAdjacent() < (_maxMerge)){ // potential 3-hit cluster, can use simple average merge.

        vector <Voxel_tpc*>* hitsToMerge = new vector <Voxel_tpc*>;

        int clusterSize = seed_hit->clusterFind(hitsToMerge);

        if( clusterSize <= _maxMerge ){ // merge cluster
          seed_hit->setIsMerged();
          writeMergedVoxelsToHit(hitsToMerge);
        }
        delete hitsToMerge;
      }
    }
  }

  int numberOfHits(0);
  // count up the number of hits merged or lost
  for (unsigned int i = 0; i<_tpcRowHits.size(); ++i){
    row_hits = _tpcRowHits.at(i);
    for (unsigned int j = 0; j<row_hits.size(); ++j){
      numberOfHits++;
      Voxel_tpc* seed_hit = row_hits[j];
      if(seed_hit->IsMerged() || seed_hit->IsClusterHit() || seed_hit->getNumberOfAdjacent() > _maxMerge ) {
        ++_NRevomedHits;
        auto mcp = (_tpcHitMap[ seed_hit ]).getMCParticle() ;
        if(mcp.isAvailable()) {
          ++_NLostPhysicsTPCHits;
          const edm4hep::Vector3f mom= mcp.getMomentum() ;
          double ptSQRD = mom[0]*mom[0]+mom[1]*mom[1] ;
          if( ptSQRD > (0.2*0.2) ) ++_NLostPhysicsAbove02GeVPtTPCHits ;
          if( ptSQRD > 1.0 )  ++_NLostPhysicsAbove1GeVPtTPCHits ;
        }
      }
    }
  }

  debug() << "the number of adjacent hits is " <<  number_of_adjacent_hits << "  _doubleHitResZ " << _doubleHitResZ << endmsg;
  debug() << "number of rec_hits = "  << _NRecTPCHits << endmsg ;
  debug() << "finished row hits " << numberOfHits << " " << numberOfhitsTreated << endmsg;

  // set the parameters to decode the type information in the collection
  // for the time being this has to be done manually
  // in the future we should provide a more convenient mechanism to
  // decode this sort of meta information

  //StringVec typeNames ;
  //IntVec typeValues ;
  //typeNames.push_back( LCIO::TPCHIT ) ;
  //typeValues.push_back( 1 ) ;
  //yzhang TODO
  //_trkhitVec->parameters().setValues("TrackerHitTypeNames" , typeNames ) ;
  //_trkhitVec->parameters().setValues("TrackerHitTypeValues" , typeValues ) ;

  //// add the collection to the event
  //evt->addCollection( _trkhitVec , _TPCTrackerHitsCol ) ;
  //evt->addCollection( _relCol , _outRelColName ) ;

  // delete voxels
  for (unsigned int i = 0; i<_tpcRowHits.size(); ++i){
    vector <Voxel_tpc *>* current_row = &_tpcRowHits.at(i);
    for (unsigned int j = 0; j<current_row->size(); ++j){
      delete current_row->at(j);
    }
  }

//#ifdef DIGIPLOTS
  //_NSimTPCHitsHisto->fill(_NSimTPCHits);
  //_NBackgroundSimTPCHitsHisto->fill(_NBackgroundSimTPCHits);
  //_NPhysicsSimTPCHitsHisto->fill(_NPhysicsSimTPCHits);
  //_NPhysicsAbove02GeVSimTPCHitsHisto->fill(_NPhysicsAbove02GeVSimTPCHits);
  //_NPhysicsAbove1GeVSimTPCHitsHisto->fill(_NPhysicsAbove1GeVSimTPCHits);
  //_NRecTPCHitsHisto->fill(_NRecTPCHits);

  //_NLostPhysicsTPCHitsHisto->fill(_NLostPhysicsTPCHits);
  //_NLostPhysicsAbove02GeVPtTPCHitsHisto->fill(_NLostPhysicsAbove02GeVPtTPCHits);
  //_NLostPhysicsAbove1GeVPtTPCHitsHisto->fill(_NLostPhysicsAbove1GeVPtTPCHits);
  //_NRevomedHitsHisto->fill(_NRevomedHits);

  //_NKeptPhysicsTPCHitsHistoPercent->fill( (float)(_NPhysicsSimTPCHits-_NLostPhysicsTPCHits) / (float)_NPhysicsSimTPCHits );
  //_NKeptPhysicsAbove02GeVPtTPCHitsHistoPercent->fill( (float)(_NPhysicsAbove02GeVSimTPCHits-_NLostPhysicsAbove02GeVPtTPCHits) / (float)_NPhysicsAbove02GeVSimTPCHits );
  //_NKeptPhysicsAbove1GeVPtTPCHitsHistoPercent->fill( (float)(_NPhysicsAbove1GeVSimTPCHits-_NLostPhysicsAbove1GeVPtTPCHits) / (float)_NPhysicsAbove1GeVSimTPCHits );
//#endif

  debug() << "_NSimTPCHits = " << _NSimTPCHits << endmsg;
  debug() << "_NBackgroundSimTPCHits = " << _NBackgroundSimTPCHits << endmsg;
  debug() << "_NPhysicsSimTPCHits = " << _NPhysicsSimTPCHits << endmsg;
  debug() << "_NPhysicsAbove02GeVSimTPCHits = " << _NPhysicsAbove02GeVSimTPCHits << endmsg;
  debug() << "_NPhysicsAbove1GeVSimTPCHits = " << _NPhysicsAbove1GeVSimTPCHits << endmsg;
  debug() << "_NRecTPCHits = " << _NRecTPCHits<< endmsg;
  debug() << "_NLostPhysicsTPCHits = " << _NLostPhysicsTPCHits << endmsg;
  debug() << "_NLostPhysicsAbove02GeVPtTPCHits = " << _NLostPhysicsAbove02GeVPtTPCHits << endmsg;
  debug() << "_NLostPhysicsAbove1GeVPtTPCHits = " << _NLostPhysicsAbove1GeVPtTPCHits << endmsg;
  debug() << "_NRevomedHits = " << _NRevomedHits << endmsg;

  _nEvt++;
  //Clear the maps and the end of the event.
  _tpcHitMap.clear();
  _tpcRowHits.clear();

  delete _cellid_encoder ;
  return StatusCode::SUCCESS;
}



StatusCode TPCDigiAlg::finalize()
{

#ifdef DIGIPLOTS
  //_TREE->commit();
  //_TREE->cd("/Histograms");
  //_TREE->ls("..");

  //_TREE->close();
  //info() << "DIGICHECKPLOTS Finished" << endmsg;
#endif

  gsl_rng_free(_random);
  info() << "TPCDigiAlg::end()  " << name()
  << " processed " << _nEvt << " events in " << _nRun << " runs "
  << endl ;

  return StatusCode::SUCCESS;
}

void TPCDigiAlg::writeVoxelToHit( Voxel_tpc* aVoxel){

  const gear::TPCParameters& gearTPC = _GEAR->getTPCParameters() ;
  const gear::PadRowLayout2D& padLayout = gearTPC.getPadLayout() ;
  const gear::Vector2D padCoord = padLayout.getPadCenter(1) ;

  Voxel_tpc* seed_hit  = aVoxel;

  //  if( seed_hit->getRowIndex() > 5 ) return ;
  debug() << "==============" << endmsg;
  //store hit variables
  edm4hep::TrackerHit trkHit;// = _trkhitVec->create();
  //now the hit pos has to be smeared

  double tpcRPhiRes = seed_hit->getRPhiRes();
  double tpcZRes = seed_hit->getZRes();

  CLHEP::Hep3Vector point(seed_hit->getX(),seed_hit->getY(),seed_hit->getZ());

  double unsmearedPhi = point.phi();

  double randrp = gsl_ran_gaussian(_random,tpcRPhiRes);
  double randz =  gsl_ran_gaussian(_random,tpcZRes);

  point.setPhi( point.phi() + randrp/ point.perp() );
  point.setZ( point.z() + randz );

  // make sure the hit is not smeared beyond the TPC Max DriftLength
  if( fabs(point.z()) > gearTPC.getMaxDriftLength() ) point.setZ( (fabs(point.z()) / point.z() ) * gearTPC.getMaxDriftLength() );
  debug() << "==============" << endmsg;
  edm4hep::Vector3d pos(point.x(),point.y(),point.z());
  trkHit.setPosition(pos);
  trkHit.setEDep(seed_hit->getEDep());
  //  trkHit->setType( 500 );

  //  int side = lcio::ILDDetID::barrel ;
  //
  //  if( pos[2] < 0.0 ) side = 1 ;
  //change to Marlin's, fucd 
  //map<Voxel_tpc*,edm4hep::SimTrackerHit>::iterator it=_tpcHitMap.find(seed_hit);
  //assert(_tpcHitMap.end() != it);

  //const int celId = it->second.getCellID();
  //_cellid_encoder->setValue( celId );
  //trkHit.setCellID( _cellid_encoder->lowWord() );

  //int side   = (*_cellid_encoder)[lcio::ILDCellID0::side];
  //int layer  = (*_cellid_encoder)[lcio::ILDCellID0::layer];
  //int module = (*_cellid_encoder)[lcio::ILDCellID0::module];
  //int sensor = (*_cellid_encoder)[lcio::ILDCellID0::sensor];

  //debug() << "Hit = " << " has celId " << celId << endmsg;
  //debug() << "side = " << side << endmsg;
  //debug() << "layerNumber = " <<  layer << endmsg;
  //debug() << "moduleNumber = " << module << endmsg;
  //debug() << "sensorNumber = " << sensor << endmsg;


  (*_cellid_encoder)[ lcio::ILDCellID0::subdet ] = lcio::ILDDetID::TPC ;
  (*_cellid_encoder)[ lcio::ILDCellID0::layer  ] = seed_hit->getRowIndex() ;
  (*_cellid_encoder)[ lcio::ILDCellID0::module ] = 0 ;

  //SJA:FIXME: for now don't use side
  //  (*_cellid_encoder)[ lcio::ILDCellID0::side   ] = side ;
  (*_cellid_encoder)[ lcio::ILDCellID0::side   ] = lcio::ILDDetID::barrel ;
  trkHit.setCellID(_cellid_encoder->lowWord());
  //_cellid_encoder->setCellID( &trkHit ) ;


  // check values for inf and nan
  if( std::isnan(unsmearedPhi) || std::isinf(unsmearedPhi) || std::isnan(tpcRPhiRes) || std::isinf(tpcRPhiRes) ) {
    std::stringstream errorMsg;
    errorMsg << "\nProcessor: TPCDigiAlg \n"
    << "unsmearedPhi = "
    <<  unsmearedPhi
    << " tpcRPhiRes = "
    <<  tpcRPhiRes
    << "\n" ;
    throw errorMsg.str();
  }
  debug() << "==============" << endmsg;
  // For no error in R
  std::array<float,TRKHITNCOVMATRIX> covMat={sin(unsmearedPhi)*sin(unsmearedPhi)*tpcRPhiRes*tpcRPhiRes,
    -cos(unsmearedPhi)*sin(unsmearedPhi)*tpcRPhiRes*tpcRPhiRes,
    cos(unsmearedPhi)*cos(unsmearedPhi)*tpcRPhiRes*tpcRPhiRes,
    0.,
    0.,
    float(tpcZRes*tpcZRes)};

  trkHit.setCovMatrix(covMat);

  if( !_tpcHitMap[seed_hit].isAvailable()){
    std::stringstream errorMsg;
    errorMsg << "\nProcessor: TPCDigiAlg \n"
    << "SimTracker Pointer is NULL throwing exception\n"
    << "\n" ;
    throw errorMsg.str();
  }
  debug() << "==============" << endmsg;
  if(pos[0]*pos[0]+pos[1]*pos[1]>0.0){
    //    push back the SimTHit for this TrackerHit

    if (_use_raw_hits_to_store_simhit_pointer) {
      trkHit.addToRawHits(_tpcHitMap[seed_hit].getObjectID());
    }

    auto rel = _relCol->create();
    rel.setRec (trkHit);
    rel.setSim (_tpcHitMap[seed_hit]);
    rel.setWeight( 1.0 );

    _trkhitVec->push_back( trkHit );
    _NRecTPCHits++;
  }

  debug() << "==============" << endmsg;
#ifdef DIGIPLOTS
//  edm4hep::SimTrackerHit* theSimHit = _tpcHitMap[seed_hit];
//  double rSimSqrd = theSimHit->getPosition()[0]*theSimHit->getPosition()[0] + theSimHit->getPosition()[1]*theSimHit->getPosition()[1];
//
//  double phiSim = atan2(theSimHit->getPosition()[1],theSimHit->getPosition()[0]);
//
//  double rPhiDiff = (point.phi() - phiSim)*sqrt(rSimSqrd);
//  double rPhiPull = ((point.phi() - phiSim)*sqrt(rSimSqrd))/(sqrt((covMat[2])/(cos(point.phi())*cos(point.phi()))));
//
//  double zDiff = point.getZ() - theSimHit->getPosition()[2];
//  double zPull = zDiff/sqrt(covMat[5]);
//
//
  //_rPhiDiffHisto->fill(rPhiDiff);
  //_rPhiPullHisto->fill(rPhiPull);
  //_phiDistHisto->fill(point.phi() - phiSim);
  //_zDiffHisto->fill(zDiff);
  //_zPullHisto->fill(zPull);

  //_zSigmaVsZHisto->fill(seed_hit->getZ(),sqrt(covMat[5]));
  //_rPhiSigmaHisto->fill(sqrt((covMat[2])/(cos(point.phi())*cos(point.phi()))));
  //_zSigmaHisto->fill(sqrt(covMat[5]));
#endif
}

void TPCDigiAlg::writeMergedVoxelsToHit( vector <Voxel_tpc*>* hitsToMerge){

  const gear::TPCParameters& gearTPC = _GEAR->getTPCParameters() ;
  const gear::PadRowLayout2D& padLayout = gearTPC.getPadLayout() ;
  const gear::Vector2D padCoord = padLayout.getPadCenter(1) ;

  edm4hep::TrackerHit trkHit;// = _trkhitVec->create();

  double sumZ = 0;
  double sumPhi = 0;
  double sumEDep = 0;
  //  double R = 0;
  double lastR = 0;

  unsigned number_of_hits_to_merge = hitsToMerge->size();


  for(unsigned int ihitCluster = 0; ihitCluster < number_of_hits_to_merge; ++ihitCluster){

    sumZ += hitsToMerge->at(ihitCluster)->getZ();
    sumPhi += hitsToMerge->at(ihitCluster)->getPhi();
    sumEDep += hitsToMerge->at(ihitCluster)->getEDep();
    hitsToMerge->at(ihitCluster)->setIsMerged();
    lastR = hitsToMerge->at(ihitCluster)->getR();

    if (_use_raw_hits_to_store_simhit_pointer) {
      trkHit.addToRawHits(_tpcHitMap[hitsToMerge->at(ihitCluster)].getObjectID());
    }

    auto rel = _relCol->create();
    rel.setRec (trkHit);
    rel.setSim (_tpcHitMap[ hitsToMerge->at(ihitCluster) ]);
    rel.setWeight( float(1.0/number_of_hits_to_merge) );

  }

  double avgZ = sumZ/(hitsToMerge->size());
  double avgPhi = sumPhi/(hitsToMerge->size());

  CLHEP::Hep3Vector* mergedPoint = new CLHEP::Hep3Vector(1.0,1.0,1.0);
  mergedPoint->setPerp(lastR);
  mergedPoint->setPhi(avgPhi);
  mergedPoint->setZ(avgZ);

  //store hit variables

  // first the hit pos has to be smeared------------------------------------------------

  //FIXME: which errors should we use for smearing the merged hits ?
  //       this might be a bit large ....
  double tpcRPhiRes = _padWidth;
  double tpcZRes = _binningZ;

  CLHEP::Hep3Vector point( mergedPoint->getX(), mergedPoint->getY(), mergedPoint->getZ()  ) ;

//  double unsmearedPhi = point.phi();

  double randrp = gsl_ran_gaussian(_random,tpcRPhiRes);
  double randz =  gsl_ran_gaussian(_random,tpcZRes);

  point.setPhi( point.phi() + randrp/ point.perp() );
  point.setZ( point.z() + randz );

  // make sure the hit is not smeared beyond the TPC Max DriftLength
  if( fabs(point.z()) > gearTPC.getMaxDriftLength() ) point.setZ( (fabs(point.z()) / point.z() ) * gearTPC.getMaxDriftLength() );

  double pos[3] = {point.x(),point.y(),point.z()};

  //---------------------------------------------------------------------------------
  trkHit.setPosition(pos);
  trkHit.setEDep(sumEDep);
  //  trkHit->setType( 500 );

  // SJA:FIXME: here you can use the value 2 but not 3 which is odd as the width of the field is 1, only 0 and 1 should be allowed?
  int side = 1 ;
  int padIndex = padLayout.getNearestPad(mergedPoint->perp(),mergedPoint->phi());
  int row = padLayout.getRowNumber(padIndex);

  if( pos[2] < 0.0 ) side = 1 ;

  (*_cellid_encoder)[ lcio::ILDCellID0::subdet ] = lcio::ILDDetID::TPC ;
  (*_cellid_encoder)[ lcio::ILDCellID0::layer  ] = row ;
  (*_cellid_encoder)[ lcio::ILDCellID0::module ] = 0 ;
  ////SJA:FIXME: for now don't use side
  ////  (*_cellid_encoder)[ lcio::ILDCellID0::side   ] = side ;
  (*_cellid_encoder)[ lcio::ILDCellID0::side   ] = 0 ;

  //yzhang FIXME?
  //_cellid_encoder->setCellID( &trkHit ) ;
  trkHit.setCellID( _cellid_encoder->lowWord() );

  double phi = mergedPoint->getPhi();

  // check values for inf and nan
  if( std::isnan(phi) || std::isinf(phi) || std::isnan(tpcRPhiRes) || std::isinf(tpcRPhiRes) ) {
    std::stringstream errorMsg;
    errorMsg << "\nProcessor: TPCDigiAlg \n"
    << "phi = "
    <<  phi
    << " tpcRPhiRes = "
    <<  tpcRPhiRes
    << "\n" ;
    throw errorMsg.str();
  }

  // For no error in R
  std::array<float,TRKHITNCOVMATRIX> covMat={sin(phi)*sin(phi)*tpcRPhiRes*tpcRPhiRes,
    -cos(phi)*sin(phi)*tpcRPhiRes*tpcRPhiRes,
    cos(phi)*cos(phi)*tpcRPhiRes*tpcRPhiRes,
    0.,
    0.,
    float(tpcZRes*tpcZRes)};

  trkHit.setCovMatrix(covMat);

  if(pos[0]*pos[0]+pos[1]*pos[1]>0.0){
    //yzhang TODO
    _trkhitVec->push_back( trkHit );
    ++_nRechits;
  } else {
    //???yzhang TODO
    //delete trkHit;
  }

  delete mergedPoint;

}


#ifdef DIGIPLOTS
//void TPCDigiAlg::plotHelixHitResidual( edm4hep::MCParticle *mcp, CLHEP::Hep3Vector *thisPoint){
//
//  const double bField = _GEAR->getBField().at( gear::Vector3D( 0., 0., 0.) ).z() ;
//  const double FCT = 2.99792458E-4;
//  double charge = mcp->getCharge();
//  const double *mom = mcp->getMomentum();
//  double pt = sqrt(mom[0]*mom[0]+mom[1]*mom[1]);
//  double radius = pt / (FCT*bField);
//  double tanLambda = mom[2]/pt;
//  double omega = charge / radius;
//
//  if(pt>1.0) {
//
//    //FIXME SJA: this is only valid for tracks from the IP and should be done correctly for non prompt tracks
//    double Z0 = 0.;
//    double D0  = 0.;
//
//    LCVector3D refPoint(0.,0.,0);
//
//    SimpleHelix* helix = new SimpleHelix(D0,
//                                         atan2(mom[1],mom[0]),
//                                         omega,
//                                         Z0,
//                                         tanLambda,
//                                         refPoint);
//
//
//    // an almost "infinite" cylinder in z
//    LCVector3D startCylinder(0.,0.,-1000000.0);
//    LCVector3D endCylinder(0.,0.,1000000.0);
//    bool endplane=true;
//
//    LCCylinder cylinder(startCylinder,endCylinder,thisPoint->perp(),endplane);
//
//    bool pointExists = false;
//
//    double pathlength = helix->getIntersectionWithCylinder( cylinder, pointExists);
//
//    LCErrorMatrix* errors = new LCErrorMatrix();
//
//    if(pointExists){
//
//      LCVector3D intersection = helix->getPosition(pathlength, errors);
//
//      double intersectionPhi = atan2(intersection[1],intersection[0]);
//      double residualRPhi = ((intersectionPhi-thisPoint->phi())) ;
//      _ResidualsRPhiHisto->fill(residualRPhi);
//
//    }
//
//    delete errors;
//    delete helix;
//
//    const gear::TPCParameters& gearTPC = _GEAR->getTPCParameters() ;
//    const gear::PadRowLayout2D& padLayout = gearTPC.getPadLayout() ;
//
//    int row = padLayout.getRowNumber(padLayout.getNearestPad(thisPoint->perp(),thisPoint->phi()));
//    int pad = padLayout.getNearestPad(thisPoint->perp(),thisPoint->phi());
//
//    double rHit_diff = thisPoint->perp()
//    - padLayout.getPlaneExtent()[0]
//    - (( row + 0.5 )
//       * padLayout.getPadHeight(pad)) ;
//
//    _radiusCheckHisto->fill(rHit_diff);
//
//    //      info() << "$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$" << endmsg;
//    //      info() << "thisPoint->perp() = " << thisPoint->perp() << endmsg;
//    //      info() << "TPC Sensitive rMin = " << padLayout.getPlaneExtent()[0] << endmsg;
//    //      info() << "Row number + 0.5 = " <<  row + 0.5 << endmsg;
//    //      info() << "Pad Height = " <<  padLayout.getPadHeight(pad) << endmsg;
//    //      info() << "Row Height = " <<   padLayout.getRowHeight(row) << endmsg;
//    //      info() << "R_hit - TPC Rmin - ((RowIndex + 0.5 )* padheight) = " << rHit_diff << endmsg;
//
//  }
//  return;
//}
#endif


double TPCDigiAlg::getPadPhi( CLHEP::Hep3Vector *thisPoint, CLHEP::Hep3Vector* firstPoint, CLHEP::Hep3Vector* middlePoint, CLHEP::Hep3Vector* lastPoint){

  CLHEP::Hep2Vector firstPointRPhi(firstPoint->x(),firstPoint->y());
  CLHEP::Hep2Vector middlePointRPhi(middlePoint->x(),middlePoint->y());
  CLHEP::Hep2Vector lastPointRPhi(lastPoint->x(),lastPoint->y());

  // check that the points are not the same, at least at the level of a tenth of a micron
  if( (fabs( firstPointRPhi.x() - middlePointRPhi.x() ) < 1.e-05  && fabs( firstPointRPhi.y() - middlePointRPhi.y() ) < 1.e-05)
     ||
     (fabs( middlePointRPhi.x() - lastPointRPhi.x() ) < 1.e-05  && fabs( middlePointRPhi.y() - lastPointRPhi.y() ) < 1.e-05)
     ||
     (fabs( firstPointRPhi.x() - lastPointRPhi.x() ) < 1.e-05  && fabs( firstPointRPhi.y() - lastPointRPhi.y() ) < 1.e-05)
     ) {

    warning() << " TPCDigiAlg::getPadPhi "
    << "2 of the 3 SimTracker hits passed to Circle Fit are the same hit taking pad phi as PI/2\n"
    << " firstPoint->x() "  << firstPoint->x()
    << " firstPoint->y() "  << firstPoint->y()
    << " firstPoint->z() "  << firstPoint->z()
    << " middlePoint->x() "  << middlePoint->x()
    << " middlePoint->y() "  << middlePoint->y()
    << " middlePoint->z() "  << middlePoint->z()
    << " lastPoint->x() "  << lastPoint->x()
    << " lastPoint->y() "  << lastPoint->y()
    << " lastPoint.z() "  << lastPoint->z()
    << std::endl ;

    return CLHEP::twopi/4.0 ;

  }



  Circle theCircle(&firstPointRPhi, &middlePointRPhi, &lastPointRPhi);

  double localPhi = atan2((thisPoint->y() - theCircle.GetCenter()->y()) ,(thisPoint->x() - theCircle.GetCenter()->x())) + (CLHEP::twopi/4.0) ;

  if(localPhi>CLHEP::twopi) localPhi=localPhi - CLHEP::twopi;
  if(localPhi<0.0) localPhi=localPhi + CLHEP::twopi;
  if(localPhi>CLHEP::pi) localPhi = localPhi - CLHEP::pi ;

  double pointPhi = thisPoint->phi();

  if(pointPhi>CLHEP::twopi) pointPhi=pointPhi - CLHEP::twopi;
  if(pointPhi<0.0) pointPhi=pointPhi + CLHEP::twopi;
  if(pointPhi>CLHEP::pi) pointPhi = pointPhi - CLHEP::pi;

  double padPhi = fabs(pointPhi - localPhi);

  // check that the value returned is reasonable
  if( std::isnan(padPhi) || std::isinf(padPhi) ) {
    std::stringstream errorMsg;
    errorMsg << "\nProcessor: TPCDigiAlg \n"
    << "padPhi = "
    <<  padPhi
    << "\n" ;
    throw errorMsg.str();
  }

  return padPhi;

}

double TPCDigiAlg::getPadTheta(CLHEP::Hep3Vector* firstPoint, CLHEP::Hep3Vector* middlePoint, CLHEP::Hep3Vector* lastPoint){

  // Calculate thetaPad for current hit
  CLHEP::Hep2Vector firstPointRPhi(firstPoint->x(),firstPoint->y()) ;
  CLHEP::Hep2Vector middlePointRPhi(middlePoint->x(),middlePoint->y());
  CLHEP::Hep2Vector lastPointRPhi(lastPoint->x(),lastPoint->y());

  // check that the points are not the same, at least at the level of a tenth of a micron
  if( (fabs( firstPointRPhi.x() - middlePointRPhi.x() ) < 1.e-05  && fabs( firstPointRPhi.y() - middlePointRPhi.y() ) < 1.e-05)
     ||
     (fabs( middlePointRPhi.x() - lastPointRPhi.x() ) < 1.e-05  && fabs( middlePointRPhi.y() - lastPointRPhi.y() ) < 1.e-05)
     ||
     (fabs( firstPointRPhi.x() - lastPointRPhi.x() ) < 1.e-05  && fabs( firstPointRPhi.y() - lastPointRPhi.y() ) < 1.e-05)
     ) {

    warning() << " TPCDigiAlg::getPadTheta "
    << "2 of the 3 SimTracker hits passed to Circle Fit are the same hit taking pad phi as PI/2\n"
    << " firstPoint->x() "  << firstPoint->x()
    << " firstPoint->y() "  << firstPoint->y()
    << " firstPoint->z() "  << firstPoint->z()
    << " middlePoint->x() "  << middlePoint->x()
    << " middlePoint->y() "  << middlePoint->y()
    << " middlePoint->z() "  << middlePoint->z()
    << " lastPoint->x() "  << lastPoint->x()
    << " lastPoint->y() "  << lastPoint->y()
    << " lastPoint.z() "  << lastPoint->z()
    << endmsg;

    return CLHEP::twopi/4.0 ;

  }


  Circle theCircle(&firstPointRPhi, &middlePointRPhi, &lastPointRPhi);

  double deltaPhi = firstPoint->deltaPhi(*lastPoint);

  double pathlength = fabs(deltaPhi) *  theCircle.GetRadius();

  double padTheta = atan ( pathlength / fabs(lastPoint->z() - firstPoint->z()) ) ;

  double pathlength1 = 2.0 * asin( ( sqrt (
                                           (middlePointRPhi.x() - firstPointRPhi.x()) * (middlePointRPhi.x()-firstPointRPhi.x())
                                           +
                                           (middlePointRPhi.y()-firstPointRPhi.y()) * (middlePointRPhi.y()-firstPointRPhi.y())
                                           ) / 2.0 ) / theCircle.GetRadius() ) * theCircle.GetRadius()  ;


  double pathlength2 = 2.0 * asin( ( sqrt (
                                           (lastPointRPhi.x()-middlePointRPhi.x()) * (lastPointRPhi.x()-middlePointRPhi.x())
                                           +
                                           (lastPointRPhi.y()-middlePointRPhi.y()) * (lastPointRPhi.y()-middlePointRPhi.y())
                                           ) / 2.0 ) / theCircle.GetRadius() ) * theCircle.GetRadius()  ;


  padTheta = atan ((fabs(pathlength1 + pathlength2)) / (fabs(lastPoint->z() - firstPoint->z())) ) ;

  // check that the value returned is reasonable
  if( std::isnan(padTheta) || std::isinf(padTheta) ) {
    std::stringstream errorMsg;
    errorMsg << "\nProcessor: TPCDigiAlg \n"
    << "padTheta = "
    <<  padTheta
    << "\n" ;
    throw errorMsg.str();
  }

  return padTheta;

}
