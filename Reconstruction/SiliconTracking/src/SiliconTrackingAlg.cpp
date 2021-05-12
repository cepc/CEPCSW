#include "SiliconTrackingAlg.h"
#include "GearSvc/IGearSvc.h"
#include "EventSeeder/IEventSeeder.h"
#include "TrackSystemSvc/ITrackSystemSvc.h"
#include "DataHelper/Navigation.h"
#include "DataHelper/TrackerHitHelper.h"
#include "edm4hep/MCParticle.h"
#include "edm4hep/TrackerHit.h"
//#include "edm4hep/TrackerHitPlane.h"
#include "edm4hep/Track.h"
#include "edm4hep/TrackState.h"

#include <iostream>
#include <algorithm>
#include <cmath>
#include <climits>

#include <gear/GEAR.h>
#include <gear/GearMgr.h>
#include <gear/GearParameters.h>
#include <gear/VXDLayerLayout.h>
#include <gear/VXDParameters.h>
#include "gear/FTDLayerLayout.h"
#include "gear/FTDParameters.h"

#include <gear/BField.h>

#include <UTIL/BitField64.h>
#include <UTIL/BitSet32.h>
#include <UTIL/ILDConf.h>

#include "TrackSystemSvc/MarlinTrkUtils.h"
#include "TrackSystemSvc/HelixTrack.h"
#include "TrackSystemSvc/HelixFit.h"
#include "TrackSystemSvc/IMarlinTrack.h"

//#include "TrackSystemSvc/MarlinTrkDiagnostics.h"
//#ifdef MARLINTRK_DIAGNOSTICS_ON
//#include "TrackSystemSvc/DiagnosticsController.h"
//#endif

//#include "MarlinCED.h"

//#include "marlin/AIDAProcessor.h"

//---- ROOT -----
#include "TH1F.h"
#include "TH2F.h"

using namespace edm4hep ;
//using namespace marlin ;
using namespace MarlinTrk ;

using std::min;
using std::max;
using std::abs;

const int SiliconTrackingAlg::_output_track_col_quality_GOOD = 1;
const int SiliconTrackingAlg::_output_track_col_quality_FAIR = 2;
const int SiliconTrackingAlg::_output_track_col_quality_POOR = 3;

const double SiliconTrackingAlg::TWOPI = 2*M_PI;

DECLARE_COMPONENT( SiliconTrackingAlg )

SiliconTrackingAlg::SiliconTrackingAlg(const std::string& name, ISvcLocator* svcLoc)
: GaudiAlgorithm(name, svcLoc) {

  //_description = "Pattern recognition in silicon trackers";
  
  _fastfitter = new MarlinTrk::HelixFit();
  
  _encoder = new UTIL::BitField64(lcio::ILDCellID0::encoder_string);
  
  _petalBasedFTDWithOverlaps = false;
  
  // zero triplet counters
  _ntriplets = _ntriplets_good = _ntriplets_2MCP = _ntriplets_3MCP = _ntriplets_1MCP_Bad = _ntriplets_bad = 0;

  // Input Collections
  // ^^^^^^^^^^^^^^^^^
  declareProperty("HeaderCol", _headerColHdl);
  declareProperty("MCParticleCollection", _inMCColHdl, "Handle of the Input MCParticle collection");
  declareProperty("VTXHitCollection", _inVTXColHdl, "Handle of the Input VTX TrackerHits collection");
  declareProperty("FTDPixelHitCollection", _inFTDPixelColHdl, "Handle of the Input FTD TrackerHits collection");
  declareProperty("FTDSpacePointCollection", _inFTDSpacePointColHdl, "Handle of the Input FTD SpacePoints collection");
  declareProperty("SITHitCollection", _inSITColHdl, "Handle of the Input SIT TrackerHits collection");
  declareProperty("SITRawHitCollection", _inSITRawColHdl, "Handle of the SIT SpacePoints raw hit collection collection");
  declareProperty("FTDRawHitCollection", _inFTDRawColHdl, "Handle of the FTD SpacePoints raw hit collection collection");
  // Output Collections
  // ^^^^^^^^^^^^^^^^^^
  declareProperty("SiTrackCollection", _outColHdl, "Handle of the SiTrack output collection");
  //declareProperty("TrkHitRelCollection", _outRelColHdl, "Handle of TrackerHit Track relation collection");
  // Steering parameters
  // ^^^^^^^^^^^^^^^^^^^
  
  _output_track_col_quality = _output_track_col_quality_GOOD;
}



StatusCode  SiliconTrackingAlg::initialize() { 
  
  _nRun = -1 ;
  _nEvt = 0 ;
  //printParameters() ;
  
  // set up the geometery needed by KalTest
  //FIXME: for now do KalTest only - make this a steering parameter to use other fitters
  auto _trackSystemSvc = service<ITrackSystemSvc>("TrackSystemSvc");
  if ( !_trackSystemSvc ) {
    error() << "Failed to find TrackSystemSvc ..." << endmsg;
    return StatusCode::FAILURE;
  }
  _trksystem =  _trackSystemSvc->getTrackSystem(this);
  
  if( _trksystem == 0 ){
    error() << "Cannot initialize MarlinTrkSystem of Type: KalTest" <<endmsg;
    return StatusCode::FAILURE;
  }
  
  _trksystem->setOption( IMarlinTrkSystem::CFG::useQMS,        _MSOn ) ;
  _trksystem->setOption( IMarlinTrkSystem::CFG::usedEdx,       _ElossOn) ;
  _trksystem->setOption( IMarlinTrkSystem::CFG::useSmoothing,  _SmoothOn) ;
  _trksystem->init() ;
  debug() << "TrackSystem builded at address=" << _trksystem << endmsg;
#ifdef MARLINTRK_DIAGNOSTICS_ON
  
  void * dcv = _trksystem->getDiagnositicsPointer();
  DiagnosticsController* dc = static_cast<DiagnosticsController*>(dcv);
  dc->init(_MarlinTrkDiagnosticsName,_MarlinTrkDiagnosticsName, _runMarlinTrkDiagnostics);
  
#endif
  
  if(setupGearGeom()==StatusCode::FAILURE) return StatusCode::FAILURE;
  
  if (_useSIT == 0)
    _nLayers = _nLayersVTX;
  else 
    _nLayers = _nLayersVTX + _nLayersSIT;
  
  // initialise the container to have separate vectors for up to _nHitsChi2 hits.
  _tracksWithNHitsContainer.resize(_nHitsChi2);
  
  _dPhi = TWOPI/_nDivisionsInPhi;
  _dTheta = 2.0/_nDivisionsInTheta;
  _dPhiFTD = TWOPI/_nDivisionsInPhiFTD;
  // I leave this for the moment, but 0.3 is c/1e9.
  // For the cut it does not make too much of a difference
  double cutOnR = _cutOnPt/(0.3*_bField);
  cutOnR = 1000.*cutOnR;
  _cutOnOmega = 1/cutOnR;
  
  _output_track_col_quality = 0;
  
  return GaudiAlgorithm::initialize();
}

StatusCode SiliconTrackingAlg::execute(){ 
  Navigation::Instance()->Initialize();
  //_current_event = evt;
  //_allHits.reserve(1000);
  auto trkCol = _outColHdl.createAndPut();

  _output_track_col_quality = _output_track_col_quality_GOOD;
  
  // zero triplet counters
  _ntriplets = _ntriplets_good = _ntriplets_2MCP = _ntriplets_3MCP = _ntriplets_1MCP_Bad = _ntriplets_bad = 0;
  
  // Clearing the working containers from the previous event
  // FIXME: partly done at the end of the event, in CleanUp. Make it consistent.
  //_tracksWithNHitsContainer.clear();
  //_trackImplVec.clear();
  
  //_colTrackerHits.clear();
  //_colNamesTrackerHits.clear();
  //auto header = _headerColHdl.get()->at(0);
  //int evtNo = header.getEventNumber();
  //int runNo = header.getRunNumber();
  //debug() << "Processing Run[" << runNo << "]::Event[" << evtNo << "]" << endmsg;

  _trackImplVec.reserve(100);

  int successVTX = InitialiseVTX();
  //int successFTD = 0;
  int successFTD = InitialiseFTD();
  if (successVTX == 1) {
    
    debug() << "      phi          theta        layer      nh o :   m :   i  :: o*m*i " << endmsg; 
    
    for (int iPhi=0; iPhi<_nDivisionsInPhi; ++iPhi) { 
      for (int iTheta=0; iTheta<_nDivisionsInTheta;++iTheta) {
        ProcessOneSector(iPhi,iTheta); // Process one VXD sector     
      }
    }
    
    debug() << "End of Processing VXD and SIT sectors" << endmsg;
    
  }
  
  if (successFTD == 1) {
    debug() << "      phi          side        layer      nh o :   m :   i  :: o*m*i " << endmsg;
    TrackingInFTD(); // Perform tracking in the FTD
    debug() << "End of Processing FTD sectors" << endmsg;
  }

  //if(0){
  if (successVTX == 1 || successFTD == 1) {
    //if (successVTX == 1 ) {
    for (int nHits = _nHitsChi2; nHits >= 3 ;// the three is hard coded, sorry.
         // It's the minimal number to form a track
         nHits--) {
      Sorting( _tracksWithNHitsContainer.getTracksWithNHitsVec( nHits ) );
      
    }
    debug() <<  "End of Sorting " << endmsg;
        
    for (int nHits = _nHitsChi2; nHits >= 3 ;// the three is hard coded, sorry.
         // It's the minimal number to form a track
         nHits--) {
      
      TrackExtendedVec &tracksWithNHits = _tracksWithNHitsContainer.getTracksWithNHitsVec( nHits );
      for (TrackExtendedVec::iterator trackIter = tracksWithNHits.begin();
           trackIter < tracksWithNHits.end(); trackIter++) {
        CreateTrack( *trackIter );
      }
      debug() <<  "End of creating "<< nHits << " hits tracks " << endmsg;
    }
    
    if (_attachFast == 0) {
      if(successVTX) AttachRemainingVTXHitsSlow();
      if(successFTD) AttachRemainingFTDHitsSlow();
    }
    else {
      if(successVTX) AttachRemainingVTXHitsFast();
      if(successFTD) AttachRemainingFTDHitsFast();
    }
    
    debug() <<  "End of picking up remaining hits " << endmsg;

    //edm4hep::TrackCollection* trkCol = nullptr; 
    //edm4hep::LCRelationCollection* relCol = nullptr;
    //auto trkCol = _outColHdl.createAndPut();
    //auto relCol = _outRelColHdl.createAndPut();
    /*
    LCCollectionVec * trkCol = new LCCollectionVec(LCIO::TRACK);
    // if we want to point back to the hits we need to set the flag
    LCFlagImpl trkFlag(0) ;
    trkFlag.setBit( LCIO::TRBIT_HITS ) ;
    trkCol->setFlag( trkFlag.getFlag()  ) ;
    
    LCCollectionVec * relCol = NULL;
    */
    //FinalRefit(trkCol, relCol);
    FinalRefit(trkCol);
    debug() << "FinalRefit finish" << endmsg;
    // set the quality of the output collection
    switch (_output_track_col_quality) {

      case _output_track_col_quality_FAIR:
        //trkCol->parameters().setValue( "QualityCode" , "Fair"  ) ;
        break;

      case _output_track_col_quality_POOR:
        //trkCol->parameters().setValue( "QualityCode" , "Poor"  ) ;
        break;
        
      default:
        //trkCol->parameters().setValue( "QualityCode" , "Good"  ) ;
        break;
    }
    /*
    if (_UseEventDisplay) {
      this->drawEvent();
    }
    */
  }
  
  // fill event based histogram
  /*
  if (_createDiagnosticsHistograms) {
  
    // triplet histos
    _histos->fill1D(DiagnosticsHistograms::hntriplets, _ntriplets);
    _histos->fill1D(DiagnosticsHistograms::hntriplets_good, _ntriplets_good);
    _histos->fill1D(DiagnosticsHistograms::hntriplets_2MCP, _ntriplets_2MCP);
    _histos->fill1D(DiagnosticsHistograms::hntriplets_3MCP, _ntriplets_3MCP);
    _histos->fill1D(DiagnosticsHistograms::hntriplets_1MCP_Bad, _ntriplets_1MCP_Bad);
    _histos->fill1D(DiagnosticsHistograms::hntriplets_bad, _ntriplets_bad);
  
  }
  */
  /*
  const edm4hep::MCParticleCollection* mcCol = nullptr;
  try{
    mcCol =_inMCColHdl.get();
  }
  catch(...){
  }
  if(mcCol){
    int id = 0;
    for(auto mcP : *mcCol){
      float pos[3];
      float mom[3];
      pos[0] = mcP.vertex()[0];
      pos[1] = mcP.vertex()[1];
      pos[2] = mcP.vertex()[2];
      mom[0] = mcP.momentum()[0];
      mom[1] = mcP.momentum()[1];
      mom[2] = mcP.momentum()[2];
      float           charge = mcP.getCharge();
      HelixClass helix;
      helix.Initialize_VP(pos,mom,charge,_bField);
      float d0 = helix.getD0();
      float z0 = helix.getZ0();
      float omega = helix.getOmega();
      float phi0 = helix.getPhi0();
      float tanLambda = helix.getTanLambda();
      std::cout <<"MCParticle: " << evtNo << " " << id << " " << sqrt(mom[0]*mom[0]+mom[1]*mom[1]) << " " << acos(mom[2]/sqrt(mom[0]*mom[0]+mom[1]*mom[1]+mom[2]*mom[2]))
		<< " " << atan2(mom[1],mom[0])
		<< " " << d0 << " " << phi0 << " " << omega << " " << z0 << " " << tanLambda << " " << mcP.vertex() << std::endl; 
      id++;
    }
  }
  
  const edm4hep::TrackCollection* trkCol = nullptr;
  try{
    trkCol = _outColHdl.get();
  }
  catch(...){
  }
  if(trkCol){
    int id = 0;
    for(auto track : *trkCol){
      int nstate = track->trackStates_size();
      for(int i=0;i<nstate;i++){
	edm4hep::TrackState trkState = track->getTrackStates(i);
	if(trkState.location != 1) continue;
	HelixClass helix_final;
        helix_final.Initialize_Canonical(trkState.phi,trkState.D0,trkState.Z0,trkState.omega,trkState.tanLambda,_bField);
        float trkPx = helix_final.getMomentum()[0];
        float trkPy = helix_final.getMomentum()[1];
	float trkPz = helix_final.getMomentum()[2];
        float trkPt = sqrt(trkPx*trkPx+trkPy*trkPy);
	std::cout << "Track parameter: " << evtNo << " " << id << " " << trkPt << " " << acos(trkPz/sqrt(trkPt*trkPt+trkPz*trkPz)) << " " << atan2(trkPy,trkPx)
		  << " " << trkState.D0 << " " << trkState.phi << " " << trkState.omega << " " << trkState.Z0 << " " << trkState.tanLambda
		  << " " << sqrt(trkState.covMatrix[0]) << " " << sqrt(trkState.covMatrix[2]) << " " << sqrt(trkState.covMatrix[5])
		  << " " << sqrt(trkState.covMatrix[9]) << " " << sqrt(trkState.covMatrix[14]) << " " << std::endl;
	id++;
	break;
      }
    }
  }
  */
  CleanUp();
  debug() << "Event is done " << endmsg;
  _nEvt++;
  return StatusCode::SUCCESS;
}


void SiliconTrackingAlg::CleanUp() {
  
  _tracksWithNHitsContainer.clear();
  
  for (int il=0;il<_nLayers;++il) {
    for (int ip=0;ip<_nDivisionsInPhi;++ip) {
      for (int it=0;it<_nDivisionsInTheta; ++it) {
        unsigned int iCode = il + _nLayers*ip + _nLayers*_nDivisionsInPhi*it;      
        
        if( iCode >= _sectors.size()){          
          error() << "iCode index out of range: iCode =   " << iCode << " _sectors.size() = " << _sectors.size() << " exit(1) called from file " << __FILE__ << " line " << __LINE__<< endmsg;
          continue;
        }
        
        TrackerHitExtendedVec& hitVec = _sectors[iCode];
        int nH = int(hitVec.size());
        for (int iH=0; iH<nH; ++iH) {
          TrackerHitExtended * hit = hitVec[iH];
          delete hit;
        }
      }
    }
  }
  
  for (int iS=0;iS<2;++iS) {
    for (unsigned int layer=0;layer<_nlayersFTD;++layer) {
      for (int ip=0;ip<_nDivisionsInPhiFTD;++ip) {
        unsigned int iCode = iS + 2*layer + 2*_nlayersFTD*ip;
        
        if( iCode >= _sectorsFTD.size()){
          error() << "iCode index out of range: iCode =   " << iCode << " _sectorsFTD.size() = " << _sectorsFTD.size() << " exit(1) called from file " << __FILE__ << " line " << __LINE__<< endmsg;
          continue;
        }
        
        TrackerHitExtendedVec& hitVec = _sectorsFTD[iCode];
        int nH = int(hitVec.size());
        for (int iH=0; iH<nH; ++iH) {
          TrackerHitExtended * hit = hitVec[iH];
          delete hit;
        }
      }
    }
  }
  _trackImplVec.clear();
}

int SiliconTrackingAlg::InitialiseFTD() {
  
  int success = 1;
  
  _nTotalFTDHits = 0;
  _sectorsFTD.clear();
  _sectorsFTD.resize(2*_nlayersFTD*_nDivisionsInPhiFTD);
  
  // Reading in FTD Pixel Hits Collection
  //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  const edm4hep::TrackerHitCollection* hitFTDPixelCol = nullptr;
  try {
    hitFTDPixelCol = _inFTDPixelColHdl.get();
  }
  catch ( GaudiException &e ) {
    debug() << "Collection " << _inFTDPixelColHdl.fullKey() << " is unavailable in event " << _nEvt << endmsg;
    //success = 0;
  }
  
  if(hitFTDPixelCol){
    //LCCollection * hitCollection = evt->getCollection(_FTDPixelHitCollection.c_str());
    
    //_colNamesTrackerHits[hitCollection] = _FTDPixelHitCollection;    
    // _colTrackerHits.push_back(hitCollection);
    
    int nelem = hitFTDPixelCol->size();
    
    debug() << "Number of FTD Pixel Hits = " << nelem << endmsg;
    _nTotalFTDHits = nelem;
    
    //for (int ielem=0; ielem<nelem; ++ielem) {
    for(auto hit : *hitFTDPixelCol){
      if ( UTIL::BitSet32( hit.getType() )[ UTIL::ILDTrkHitTypeBit::ONE_DIMENSIONAL ] ) continue;
      //dm4hep::ConstTrackerHit hit = hitFTDPixelCol->at(ielem);
      TrackerHitExtended * hitExt = new TrackerHitExtended( hit );
      //gear::Vector3D U(1.0,hit->getU()[1],hit->getU()[0],gear::Vector3D::spherical);
      //gear::Vector3D V(1.0,hit->getV()[1],hit->getV()[0],gear::Vector3D::spherical);
      gear::Vector3D U(1.0,hit.getCovMatrix()[1],hit.getCovMatrix()[0],gear::Vector3D::spherical);
      gear::Vector3D V(1.0,hit.getCovMatrix()[4],hit.getCovMatrix()[3],gear::Vector3D::spherical);
      gear::Vector3D Z(0.0,0.0,1.0);
      
      const float eps = 1.0e-07;
      // V must be the global z axis 
      if( fabs(V.dot(Z)) > eps ) {
	error() << "SiliconTrackingAlg: FTD Hit measurment vectors V is not in the global X-Y plane. \n\n  exit(1) called from file " << __FILE__ << " and line " << __LINE__ << endmsg;
	exit(1);
      }
      
      if( fabs(U.dot(Z)) > eps ) {
	error() << "SiliconTrackingAlg: FTD Hit measurment vectors U is not in the global X-Y plane. \n\n exit(1) called from file " << __FILE__ << " and line " << __LINE__ << endmsg;
	exit(1);
      }
      // SJA:FIXME Here dU and dV are almost certainly dX and dY ... should test ...
      //double point_res_rphi = sqrt( hit->getdU()*hit->getdU() + hit->getdV()*hit->getdV() );
      double point_res_rphi = sqrt( hit.getCovMatrix()[2]*hit.getCovMatrix()[2] + hit.getCovMatrix()[5]*hit.getCovMatrix()[5] );
      hitExt->setResolutionRPhi( point_res_rphi );
      
      // SJA:FIXME why is this needed? 
      hitExt->setResolutionZ(0.1);
      
      // type is now only used in one place where it is set to 0 to reject hits from a fit, set to INT_MAX to try and catch any missuse
      hitExt->setType(int(INT_MAX));
      // det is no longer used set to INT_MAX to try and catch any missuse
      hitExt->setDet(int(INT_MAX));
    
      double pos[3];
      
      for (int i=0; i<3; ++i) {
	pos[i] = hit.getPosition()[i];
      }
      
      double Phi = atan2(pos[1],pos[0]);
      if (Phi < 0.) Phi = Phi + TWOPI;
      // get the layer number
      unsigned int layer = static_cast<unsigned int>(getLayerID(hit));
      unsigned int petalIndex = static_cast<unsigned int>(getModuleID(hit));
      
      if ( _petalBasedFTDWithOverlaps == true ) {
	
	// as we are dealing with staggered petals we will use 2*nlayers in each directions +/- z
	// the layers will follow the even odd numbering of the petals 
	if ( petalIndex % 2 == 0 ) {
	  layer = 2*layer;
	}
	else {
	  layer = 2*layer + 1;
	}
	//debug() << "_petalBasedFTDWithOverlaps = true layer->2*layer,2layer+1 according to petalIndex = " << petalIndex << endmsg;
      }
      if (layer >= _nlayersFTD) {
	fatal() << "SiliconTrackingAlg => fatal error in FTD : layer is outside allowed range : " << layer << " number of layers = " << _nlayersFTD <<  endmsg;
	exit(1);
      }
      
      int iPhi = int(Phi/_dPhiFTD);
      
      int side = getSideID(hit);
      int iSemiSphere = 0;
    
      if (side > 0) 
	iSemiSphere = 1;
      
      int iCode = iSemiSphere + 2*layer + 2*_nlayersFTD*iPhi;
      _sectorsFTD[iCode].push_back( hitExt );
      
      debug() << " FTD Pixel Hit added : @ " << pos[0] << " " << pos[1] << " " << pos[2] << " drphi " << hitExt->getResolutionRPhi() << " dz " << hitExt->getResolutionZ() << "  iPhi = " << iPhi << " iSemiSphere "  << iSemiSphere << " iCode = " << iCode << " layer = " << layer << " cellID = " << hit.getCellID() << endmsg;  
    }
  }
  // Reading out FTD SpacePoint Collection
  //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  const edm4hep::TrackerHitCollection* hitFTDSpacePointCol = nullptr;
  try {
    hitFTDSpacePointCol = _inFTDSpacePointColHdl.get();
  }
  catch ( GaudiException &e ) {
    debug() << "Collection " << _inFTDSpacePointColHdl.fullKey() << " is unavailable in event " << _nEvt << endmsg;
    //success = 0;
  }

  const edm4hep::TrackerHitCollection* rawHitCol = nullptr;
  if(hitFTDSpacePointCol){
    try{
      rawHitCol = _inFTDRawColHdl.get();
    }
    catch ( GaudiException &e ) {
      fatal() << "Collection " << _inFTDRawColHdl.fullKey() << " is unavailable in event " << _nEvt << endmsg;
      success = 0;
    }
  }
  
  if(hitFTDSpacePointCol&&rawHitCol){
    Navigation::Instance()->AddTrackerHitCollection(rawHitCol);
    //LCCollection * hitCollection = evt->getCollection(_FTDSpacePointCollection.c_str());
    
    //_colNamesTrackerHits[hitCollection] = _FTDSpacePointCollection;
    //_colTrackerHits.push_back(hitCollection);
    
    int nelem = hitFTDSpacePointCol->size();
    
    debug() << "Number of FTD SpacePoints = " << nelem << endmsg;
    _nTotalFTDHits += nelem;
    
    //for (int ielem=0; ielem<nelem; ++ielem) {
    for(auto hit : *hitFTDSpacePointCol){
    //edm4hep::ConstTrackerHit hit =  hitFTDSpacePointCol->at(ielem);
      
      TrackerHitExtended * hitExt = new TrackerHitExtended(hit);
      
      // SJA:FIXME: fudge for now by a factor of two and ignore covariance
      double point_res_rphi = 2 * sqrt( hit.getCovMatrix()[0] + hit.getCovMatrix()[2] );
      
      hitExt->setResolutionRPhi( point_res_rphi );
      
      // SJA:FIXME why is this needed? 
      hitExt->setResolutionZ(0.1);
      
      // type is now only used in one place where it is set to 0 to reject hits from a fit, set to INT_MAX to try and catch any missuse
      hitExt->setType(int(INT_MAX));
      // det is no longer used set to INT_MAX to try and catch any missuse
      hitExt->setDet(int(INT_MAX));
      
      double pos[3];
      
      for (int i=0; i<3; ++i) {
        pos[i] = hit.getPosition()[i];
      }
      
      double Phi = atan2(pos[1],pos[0]);
      if (Phi < 0.) Phi = Phi + TWOPI;
      
      // get the layer number
      unsigned int layer = static_cast<unsigned int>(getLayerID(hit));
      unsigned int petalIndex = static_cast<unsigned int>(getModuleID(hit));
      
      if ( _petalBasedFTDWithOverlaps == true ) {
        // as we are dealing with staggered petals we will use 2*nlayers in each directions +/- z
        // the layers will follow the even odd numbering of the petals 
        if ( petalIndex % 2 == 0 ) {
          layer = 2*layer;
        }
        else {
          layer = 2*layer + 1;
        }
        //debug() << "_petalBasedFTDWithOverlaps = true layer->2*layer,2layer+1 according to petalIndex = " << petalIndex << endmsg;
      }
      
      if (layer >= _nlayersFTD) {
        error() << "SiliconTrackingAlg => fatal error in FTD : layer is outside allowed range : " << layer << " number of layers = " << _nlayersFTD <<  endmsg;
        exit(1);
      }
      
      int iPhi = int(Phi/_dPhiFTD);
      
      int side = getSideID(hit);
      int iSemiSphere = 0;
      
      if (side > 0) 
        iSemiSphere = 1;
      
      int iCode = iSemiSphere + 2*layer + 2*_nlayersFTD*iPhi;
      _sectorsFTD[iCode].push_back( hitExt );
      
      debug() << " FTD SpacePoint Hit added : @ " << pos[0] << " " << pos[1] << " " << pos[2] << " drphi " << hitExt->getResolutionRPhi() << " dz " << hitExt->getResolutionZ() << "  iPhi = " << iPhi << " iSemiSphere "  << iSemiSphere << " iCode = " << iCode << " layer = " << layer << " cellID = " << hit.getCellID() << endmsg;  
      
    }
  }
  
  for (unsigned i=0; i<_sectorsFTD.size(); ++i) {
    int nhits = _sectorsFTD[i].size();
    if( nhits != 0 ) debug() << " Number of Hits in FTD Sector " << i << " = " << _sectorsFTD[i].size() << endmsg;
    if (nhits > _max_hits_per_sector) {
      for (unsigned ihit=0; ihit<_sectorsFTD[i].size(); ++ihit) {
        delete _sectorsFTD[i][ihit];
      } 
      _sectorsFTD[i].clear();
      if( nhits != 0 ) error()  << " \n ### Number of Hits in FTD Sector " << i << " = " << nhits << " : Limit is set to " << _max_hits_per_sector << " : This sector will be dropped from track search, and QualityCode set to \"Poor\" " << endmsg;
      
      _output_track_col_quality = _output_track_col_quality_POOR;
      
    }
    
  }
  if(hitFTDPixelCol==nullptr&&hitFTDSpacePointCol==nullptr) success = 0;

  debug() << "FTD initialized" << endmsg;
  return success;
}

int SiliconTrackingAlg::InitialiseVTX() {
  _nTotalVTXHits = 0;
  _nTotalSITHits = 0;
  _sectors.clear();
  _sectors.resize(_nLayers*_nDivisionsInPhi*_nDivisionsInTheta);
  int success = 1;
  // Reading out VTX Hits Collection
  //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^  
  const edm4hep::TrackerHitCollection* hitVTXCol = nullptr;
  try {
    hitVTXCol = _inVTXColHdl.get();
  }
  catch ( GaudiException &e ) {
    debug() << "Collection " << _inVTXColHdl.fullKey() << " is unavailable in event " << _nEvt << endmsg;
    success = 0;
  }
  if(hitVTXCol){
    //LCCollection * hitCollection = evt->getCollection(_VTXHitCollection.c_str());
    
    //_colNamesTrackerHits[hitCollection] = _VTXHitCollection;
    //_colTrackerHits.push_back(hitCollection);
    int nelem = hitVTXCol->size();
    debug() << "Number of VTX hits = " << nelem << endmsg;
    _nTotalVTXHits = nelem;
    
    for (int ielem=0; ielem<nelem; ++ielem) {
    //for(auto hit : *hitVTXCol){ 
      edm4hep::TrackerHit hit = hitVTXCol->at(ielem);
      //gear::Vector3D U(1.0,hit->getU()[1],hit->getU()[0],gear::Vector3D::spherical);
      //gear::Vector3D V(1.0,hit->getV()[1],hit->getV()[0],gear::Vector3D::spherical);
      gear::Vector3D U(1.0,hit.getCovMatrix()[1],hit.getCovMatrix()[0],gear::Vector3D::spherical);
      gear::Vector3D V(1.0,hit.getCovMatrix()[4],hit.getCovMatrix()[3],gear::Vector3D::spherical);
      gear::Vector3D Z(0.0,0.0,1.0);
      //debug() << "covMatrix : " << hit->getCovMatrix()[0] << " " << hit->getCovMatrix()[1] << endmsg;
      const float eps = 1.0e-07;
      // V must be the global z axis 
      if( fabs(1.0 - V.dot(Z)) > eps ) {
        error() << "SiliconTrackingAlg: VXD Hit measurment vectors V is not equal to the global Z axis. \n\n  exit(1) called from file " << __FILE__ << " and line " << __LINE__ << endmsg;
        exit(1);
      }
      
      if( fabs(U.dot(Z)) > eps ) {
        error() << "SiliconTrackingAlg: VXD Hit measurment vectors U is not in the global X-Y plane. \n\n exit(1) called from file " << __FILE__ << " and line " << __LINE__ << endmsg;
        exit(1);
      }
      TrackerHitExtended * hitExt = new TrackerHitExtended(hit);
      //debug() << "Saved TrackerHit id in TrackerHitExtended " << ielem << ": " << hitExt->getTrackerHit().id() << std::endl;
            
      // SJA:FIXME: just use planar res for now
      hitExt->setResolutionRPhi(hit.getCovMatrix()[2]);
      hitExt->setResolutionZ(hit.getCovMatrix()[5]);
      
      // set type is now only used in one place where it is set to 0 to reject hits from a fit, set to INT_MAX to try and catch any missuse
      hitExt->setType(int(INT_MAX));
      // det is no longer used set to INT_MAX to try and catch any missuse
      hitExt->setDet(int(INT_MAX));
      
      double pos[3];
      double radius = 0;
      
      for (int i=0; i<3; ++i) {
        pos[i] = hit.getPosition()[i];
        radius += pos[i]*pos[i];
      }
      
      radius = sqrt(radius);
      
      double cosTheta = pos[2]/radius;
      double Phi = atan2(pos[1],pos[0]);
      
      if (Phi < 0.) Phi = Phi + TWOPI;
      
      // get the layer number
      int layer = getLayerID(hit);
      
      if (layer < 0 || layer >= _nLayers) {
        error() << "SiliconTrackingAlg => fatal error in VTX : layer is outside allowed range : " << layer << endmsg;
        exit(1);
      }
      
      int iPhi = int(Phi/_dPhi);
      int iTheta = int ((cosTheta + double(1.0))/_dTheta);
      int iCode = layer + _nLayers*iPhi + _nLayers*_nDivisionsInPhi*iTheta;      
      _sectors[iCode].push_back( hitExt );
      
      debug() << " VXD Hit " <<  hit.id() << " added : @ " << pos[0] << " " << pos[1] << " " << pos[2] << " drphi " << hitExt->getResolutionRPhi() << " dz " << hitExt->getResolutionZ() << "  iPhi = " << iPhi <<  " iTheta "  << iTheta << " iCode = " << iCode << "  layer = " << layer << endmsg;  
      
    }
  }
  
  if (_useSIT > 0 ) {
    const edm4hep::TrackerHitCollection* hitSITCol = nullptr;
    try {
      hitSITCol = _inSITColHdl.get();
    }
    catch ( GaudiException &e ) {
      debug() << "Collection " << _inSITColHdl.fullKey() << " is unavailable in event " << _nEvt << endmsg;
      success = 0;
    }

    const edm4hep::TrackerHitCollection* rawHitCol = nullptr;
    if(hitSITCol){
      try{
	rawHitCol = _inSITRawColHdl.get();
      }
      catch ( GaudiException &e ) {
	warning() << "Collection " << _inSITRawColHdl.fullKey() << " is unavailable in event " << _nEvt << ", if SIT is Space Point, it needed " << endmsg;
	//success = 0;
      }
    }
    
    if(hitSITCol){
      if(rawHitCol) Navigation::Instance()->AddTrackerHitCollection(rawHitCol);
      //debug() << "SITHitCollection pointer = " << hitSITCol << endmsg;
      int nelem = hitSITCol->size();
      
      debug() << "Number of SIT hits = " << nelem << endmsg;
      _nTotalSITHits = nelem;
      
      //ConstTrackerHit          trkhit   = 0;
      //TrackerHitPlane*     trkhit_P = 0;
      //TrackerHitZCylinder* trkhit_C = 0;
      
      double drphi(NAN);
      double dz(NAN);
      
      //for (int ielem=0; ielem<nelem; ++ielem) {
      for(auto trkhit : *hitSITCol){
        // hit could be of the following type
        // 1) TrackerHit, either ILDTrkHitTypeBit::COMPOSITE_SPACEPOINT or just standard TrackerHit
        // 2) TrackerHitPlane, either 1D or 2D
        // 3) TrackerHitZCylinder, if coming from a simple cylinder design as in the LOI
        
        // Establish which of these it is in the following order of likelyhood
        //    i)   ILDTrkHitTypeBit::ONE_DIMENSIONAL (TrackerHitPlane) Should Never Happen: SpacePoints Must be Used Instead
        //    ii)  ILDTrkHitTypeBit::COMPOSITE_SPACEPOINT (TrackerHit)
        //    iii) TrackerHitPlane (Two dimentional)
        //    iv)  TrackerHitZCylinder 
        //    v)   Must be standard TrackerHit
        
	//const edm4hep::ConstTrackerHit trkhit = hitSITCol->at(ielem);
        int layer = getLayerID(trkhit);
        
        // VXD and SIT are treated as one system so SIT layers start from _nLayersVTX
        layer = layer + _nLayersVTX;
        
        if (layer < 0 || layer >= _nLayers) {
          error() << "SiliconTrackingAlg => fatal error in SIT : layer is outside allowed range : " << layer << endmsg;
          exit(1);
        }

        // first check that we have not been given 1D hits by mistake, as they won't work here
        if ( UTIL::BitSet32( trkhit.getType() )[ UTIL::ILDTrkHitTypeBit::ONE_DIMENSIONAL ] ) {
	  fatal() << "SiliconTrackingAlg: SIT Hit cannot be of type UTIL::ILDTrkHitTypeBit::ONE_DIMENSIONAL COMPOSITE SPACEPOINTS must be use instead. \n\n  exit(1) called from file " << __FILE__ << " and line " << __LINE__ << endmsg;
	  exit(1);
        } 
        // most likely case: COMPOSITE_SPACEPOINT hits formed from stereo strip hits
        else if ( UTIL::BitSet32( trkhit.getType() )[ UTIL::ILDTrkHitTypeBit::COMPOSITE_SPACEPOINT ] ) {
          
          // SJA:FIXME: fudge for now by a factor of two and ignore covariance
          drphi =  2 * sqrt(trkhit.getCovMatrix()[0] + trkhit.getCovMatrix()[2]);         
          dz    =      sqrt(trkhit.getCovMatrix()[5]);         
          
        } 
        // or a PIXEL based SIT, using 2D TrackerHitPlane like the VXD above
	// by fucd
        //else if ( ( trkhit_P = dynamic_cast<TrackerHitPlane*>( hitCollection->getElementAt( ielem ) ) ) )  {
	else if( UTIL::BitSet32( trkhit.getType() )[ 3 ]){
          // first we need to check if the measurement vectors are aligned with the global coordinates 
          //gear::Vector3D U(1.0,trkhit_P->getU()[1],trkhit_P->getU()[0],gear::Vector3D::spherical);
          //gear::Vector3D V(1.0,trkhit_P->getV()[1],trkhit_P->getV()[0],gear::Vector3D::spherical);
	  gear::Vector3D U(1.0,trkhit.getCovMatrix()[1],trkhit.getCovMatrix()[0],gear::Vector3D::spherical);
	  gear::Vector3D V(1.0,trkhit.getCovMatrix()[4],trkhit.getCovMatrix()[3],gear::Vector3D::spherical);
          gear::Vector3D Z(0.0,0.0,1.0);
          
          const float eps = 1.0e-07;
          // V must be the global z axis 
          if( fabs(1.0 - V.dot(Z)) > eps ) {
            error() << "SiliconTrackingAlg: PIXEL SIT Hit measurment vectors V is not equal to the global Z axis. \n\n  exit(1) called from file " << __FILE__ << " and line " << __LINE__ << endmsg;
            exit(1);
          }
          
          // U must be normal to the global z axis
          if( fabs(U.dot(Z)) > eps ) {
            error() << "SiliconTrackingAlg: PIXEL SIT Hit measurment vectors U is not in the global X-Y plane. \n\n exit(1) called from file " << __FILE__ << " and line " << __LINE__ << endmsg;
            exit(1);
          }
          
          //drphi = trkhit_P->getdU();
          //dz    = trkhit_P->getdV();                                                 
          drphi = trkhit.getCovMatrix()[2];
          dz    = trkhit.getCovMatrix()[5];
        }
	
        // or a simple cylindrical design, as used in the LOI      
	/* by fucd
        else if ( ( trkhit_C = dynamic_cast<TrackerHitZCylinder*>( hitCollection->getElementAt( ielem ) ) ) ) {
          
          drphi = trkhit_C->getdRPhi();
          dz    = trkhit_C->getdZ();
          
        } 
	*/
        // this would be very unlikely, but who knows ... just an ordinary TrackerHit, which is not a COMPOSITE_SPACEPOINT
        else {
          
          // SJA:FIXME: fudge for now by a factor of two and ignore covariance
          drphi =  2 * sqrt(trkhit.getCovMatrix()[0] + trkhit.getCovMatrix()[2]);         
          dz =     sqrt(trkhit.getCovMatrix()[5]);             
          
        }
        // now that the hit type has been established carry on and create a 
        
        TrackerHitExtended * hitExt = new TrackerHitExtended(trkhit);
        
        // SJA:FIXME: just use planar res for now
        hitExt->setResolutionRPhi(drphi);
        hitExt->setResolutionZ(dz);
        
        // set type is now only used in one place where it is set to 0 to reject hits from a fit, set to INT_MAX to try and catch any missuse
        hitExt->setType(int(INT_MAX));
        // det is no longer used set to INT_MAX to try and catch any missuse
        hitExt->setDet(int(INT_MAX));
        
        double pos[3];
        double radius = 0;
        
        for (int i=0; i<3; ++i) {
          pos[i] = trkhit.getPosition()[i];
          radius += pos[i]*pos[i];
        }
        
        radius = sqrt(radius);
        
        double cosTheta = pos[2]/radius;
        double Phi = atan2(pos[1],pos[0]);
        
        if (Phi < 0.) Phi = Phi + TWOPI;
        
        int iPhi = int(Phi/_dPhi);
        int iTheta = int ((cosTheta + double(1.0))/_dTheta);
        int iCode = layer + _nLayers*iPhi + _nLayers*_nDivisionsInPhi*iTheta;      
        _sectors[iCode].push_back( hitExt );
        
        debug() << " SIT Hit " <<  trkhit.id() << " added : @ " << pos[0] << " " << pos[1] << " " << pos[2] << " drphi " << hitExt->getResolutionRPhi() << " dz " << hitExt->getResolutionZ() << "  iPhi = " << iPhi <<  " iTheta "  << iTheta << " iCode = " << iCode << "  layer = " << layer << endmsg;  
        
      }
    }    
  }
  
  
  for (unsigned i=0; i<_sectors.size(); ++i) {
    int nhits = _sectors[i].size();
    if( nhits != 0 ) debug() << " Number of Hits in VXD/SIT Sector " << i << " = " << _sectors[i].size() << endmsg;
    if (nhits > _max_hits_per_sector) {
      for (unsigned ihit=0; ihit<_sectors[i].size(); ++ihit) {
        delete _sectors[i][ihit];
      }
      _sectors[i].clear();
      if( nhits != 0 ) error()  << " \n ### Number of Hits in VXD/SIT Sector " << i << " = " << nhits << " : Limit is set to " << _max_hits_per_sector << " : This sector will be dropped from track search, and QualityCode set to \"Poor\" " << endmsg;
      
      _output_track_col_quality = _output_track_col_quality_POOR;
    }
  }
  debug() << "VXD initialized" << endmsg;
  return success; 
}

StatusCode  SiliconTrackingAlg::finalize(){
  
  delete _fastfitter ; _fastfitter = 0;
  delete _encoder ; _encoder = 0;
  //delete _trksystem ; _trksystem = 0;
  //delete _histos ; _histos = 0;
  info() << "Processed " << _nEvt << " events " << endmsg;
  info() << lcio::ILDCellID0::encoder_string << " " << UTIL::ILDCellID0::encoder_string << endmsg;
  return GaudiAlgorithm::finalize();
}


void SiliconTrackingAlg::ProcessOneSector(int iPhi, int iTheta) {
  
  int counter = 0 ;
  
  int iPhi_Up    = iPhi + 1;
  int iPhi_Low   = iPhi - 1;
  int iTheta_Up  = iTheta + 1; 
  int iTheta_Low = iTheta - 1;
  if (iTheta_Low < 0) iTheta_Low = 0;
  if (iTheta_Up  >= _nDivisionsInTheta) iTheta_Up = _nDivisionsInTheta-1;
  
  int nComb = int( _Combinations.size() / 3 ); // number of triplet combinations
                                               //  std::cout << iPhi << " " << iTheta << " " << _nEvt << endmsg;
  int iNC = 0;
  
  for (int iComb=0; iComb < nComb; ++iComb) { // loop over triplets
    
    int nLR[3];
    
    for (int iS=0; iS<3; ++iS) {
      nLR[iS] = _Combinations[iNC];
      iNC++;
    }    
    
    //std::cout << iPhi << " " << iTheta << " " << nLR[0] << " " << nLR[1] << " " << nLR[2] << " " << std::endl;
    
    // index of theta-phi bin of outer most layer
    int iCode = nLR[0] + _nLayers*iPhi +  _nLayers*_nDivisionsInPhi*iTheta;
    
    //std::cout << "size of vector = " << _sectors.size() << " iCode = " << iCode << std::endl;
    
    // get the all the hits in the outer most theta-phi bin 
    
    TrackerHitExtendedVec& hitVecOuter =  _sectors[iCode];
    
    int nHitsOuter = int(hitVecOuter.size());
    if (nHitsOuter > 0) {
      
      //std::cout << " " << iPhi << " " << iTheta << " " << nLR[0] << " " << nLR[1] << " " << nLR[2] << " size of vector = " << hitVecOuter.size() << std::endl;
      
      for (int ipMiddle=iPhi_Low; ipMiddle<iPhi_Up+1;ipMiddle++) { // loop over phi in the Middle
        
        for (int itMiddle=iTheta_Low; itMiddle<iTheta_Up+1;itMiddle++) { // loop over theta in the Middle 
          
          int iPhiMiddle = ipMiddle;
          
          // catch wrap-around
          if (ipMiddle < 0) iPhiMiddle = _nDivisionsInPhi-1;          
          if (ipMiddle >= _nDivisionsInPhi) iPhiMiddle = ipMiddle - _nDivisionsInPhi;
          
          // index of current theta-phi bin of middle layer
          iCode = nLR[1] + _nLayers*iPhiMiddle +  _nLayers*_nDivisionsInPhi*itMiddle;
          
          // get the all the hits in the current middle theta-phi bin 
          TrackerHitExtendedVec& hitVecMiddle = _sectors[iCode];
          
          int nHitsMiddle = int(hitVecMiddle.size());
          
          // determine which inner theta-phi bins to look in
          
          int iPhiLowInner = iPhi_Low;
          int iPhiUpInner = iPhi_Up;
          int iThetaLowInner = iTheta_Low;
          int iThetaUpInner = iTheta_Up;        
          
          // test to see if this is the core bin of the current search
          // if so, look into the neigboring bins in the inner layer
          if (ipMiddle == iPhi && itMiddle==iTheta) { 
            iPhiLowInner = iPhi_Low;
            iPhiUpInner  = iPhi_Up;
            iThetaLowInner = iTheta_Low;
            iThetaUpInner = iTheta_Up;
          }
          else { 
            int difP = abs(ipMiddle-iPhi); //  number of phi bins from core: can only be 1 or 0 due to hard coded 1 above
            int difT = abs(itMiddle-iTheta);// number of theta bins from core: can only be 1 or 0 due to hard coded 1 above
            int minP = min(ipMiddle,iPhi);   // min phi: core bin or current phi bin middle
            int minT = min(itMiddle,iTheta); // min theta: core bin or current theta bin middle
            int maxP = max(ipMiddle,iPhi);   // max phi: core bin or current phi bin middle
            int maxT = max(itMiddle,iTheta); // max theta: core bin or current theta bin middle
            
            if (difP==1 && difT==1) { // if the diffence is a single bin in both phi and theta : only look in the bin adjacent to the core bin  
              iPhiLowInner = minP;
              iPhiUpInner = maxP;
              iThetaLowInner = minT;
              iThetaUpInner = maxT;
            }
            if (difP==0) { // must be +/-1 theta : only look in bins adjacent to the middle bin
              iPhiLowInner = iPhi_Low;
              iPhiUpInner  = iPhi_Up;
              iThetaLowInner = minT;
              iThetaUpInner = maxT;
            }
            if (difT==0) { // must be +/-1 phi : only look in bins adjacent to the middle bin
              iPhiLowInner = minP;
              iPhiUpInner  = maxP;
              iThetaLowInner = iTheta_Low;
              iThetaUpInner = iTheta_Up;            
            }
          }               
          if (nHitsMiddle > 0) { // look into inner bins
            
            for (int ipInner=iPhiLowInner; ipInner<iPhiUpInner+1;ipInner++) { // loop over phi in the Inner
              
              for (int itInner=iThetaLowInner; itInner<iThetaUpInner+1;itInner++) { // loop over theta in the Inner 
                
                int iPhiInner = ipInner;
                
                // catch wrap-around
                if (ipInner < 0) iPhiInner = _nDivisionsInPhi-1;
                if (ipInner >= _nDivisionsInPhi) iPhiInner = ipInner - _nDivisionsInPhi;
                
                iCode = nLR[2] + _nLayers*iPhiInner +  _nLayers*_nDivisionsInPhi*itInner;
                
                // get hit for inner bin
                TrackerHitExtendedVec& hitVecInner = _sectors[iCode];
                
                int nHitsInner = int(hitVecInner.size());
                
                if (nHitsInner > 0) {
                  
                  debug() << " " 
                  << std::setw(3) << iPhi       << " "   << std::setw(3) << ipMiddle << " "      << std::setw(3) << ipInner << "   " 
                  << std::setw(3) << iTheta     << " "   << std::setw(3) << itMiddle << " "      << std::setw(3) << itInner << "  " 
                  << std::setw(3) << nLR[0]     << " "   << std::setw(3) << nLR[1]   << " "      << std::setw(3) << nLR[2]  << "     " 
                  << std::setw(3) << nHitsOuter << " : " << std::setw(3) << nHitsMiddle << " : " << std::setw(3) << nHitsInner << "  :: " 
                  << std::setw(3) << nHitsOuter*nHitsMiddle* nHitsInner << endmsg;
                  
                  // test all triplets 
                  
                  for (int iOuter=0; iOuter<nHitsOuter; ++iOuter) { // loop over hits in the outer sector
                    TrackerHitExtended * outerHit = hitVecOuter[iOuter];
                    for (int iMiddle=0;iMiddle<nHitsMiddle;iMiddle++) { // loop over hits in the middle sector
                      TrackerHitExtended * middleHit = hitVecMiddle[iMiddle];
                      for (int iInner=0;iInner<nHitsInner;iInner++) { // loop over hits in the inner sector
                        TrackerHitExtended * innerHit = hitVecInner[iInner];
                        HelixClass helix;
			// test fit to triplet
                        TrackExtended * trackAR = TestTriplet(outerHit,middleHit,innerHit,helix);
			if ( trackAR != NULL ) {
                          int nHits = BuildTrack(outerHit,middleHit,innerHit,helix,nLR[2],
                                                 iPhiLowInner,iPhiUpInner,
                                                 iThetaLowInner,iThetaUpInner,trackAR);
                          
                          _tracksWithNHitsContainer.getTracksWithNHitsVec(nHits).push_back(trackAR);
                          
                          counter ++ ;
                        }
                      } // endloop over hits in the inner sector
                    } // endloop over hits in the middle sector
                  } // endloop over hits in the outer sector
                } // endif nHitsInner > 0
              } // endloop over theta in the Inner
            } // endloop over phi in the Inner      
          } // endif nHitsMiddle > 0
        } // endloop over theta in the Middle
      } // endloop over phi in the Middle
    } // endif nHitsOuter > 0
  } // endloop over triplets
    
  //debug() << " process one sectector theta,phi " << iTheta << ", " << iPhi << "  number of loops : " << counter << endmsg  ;
}

TrackExtended * SiliconTrackingAlg::TestTriplet(TrackerHitExtended * outerHit, 
                                                       TrackerHitExtended * middleHit,
                                                       TrackerHitExtended * innerHit,
                                                       HelixClass & helix) {
  /*
   Methods checks if the triplet of hits satisfies helix hypothesis
   */
  // get the tracks already associated with the triplet
  TrackExtendedVec& trackOuterVec  = outerHit->getTrackExtendedVec();
  TrackExtendedVec& trackMiddleVec = middleHit->getTrackExtendedVec();
  TrackExtendedVec& trackInnerVec  = innerHit->getTrackExtendedVec();
  
  // check if all the hits are already assigned to a track 
  if ( (!trackOuterVec.empty())  && (!trackMiddleVec.empty()) && (!trackInnerVec.empty())) {
    
    TrackExtendedVec::const_iterator middleEndIter = trackMiddleVec.end();
    TrackExtendedVec::const_iterator outerEndIter  = trackOuterVec.end();
    TrackExtendedVec::const_iterator innerEndIter  = trackInnerVec.end();
    TrackExtendedVec::const_iterator outerBeginIter  = trackOuterVec.begin();
    TrackExtendedVec::const_iterator innerBeginIter  = trackInnerVec.begin();
    
    // loop over the tracks from the middle hit
    for (TrackExtendedVec::const_iterator middleIter = trackMiddleVec.begin();
         middleIter < middleEndIter;
         ++middleIter) {
      
      // loop over the track from the outer hit
      for (TrackExtendedVec::const_iterator outerIter = outerBeginIter;
           outerIter < outerEndIter;
           ++outerIter) {
        
        // if track from the outer and middle are not the same progress  
        if ( *outerIter != *middleIter ) continue;
        
        // loop over the tracks from the inner hit
        for (TrackExtendedVec::const_iterator innerIter = innerBeginIter;
             innerIter < innerEndIter;
             ++innerIter) {
          
          // no need to check against middle, it is idendical to outer here
          if ( *outerIter == *innerIter ) {
            // an existing track already contains all three hits
            // return a null pointer
            debug() << " TestTriplet: track " << *outerIter << " already contains all three hits: Do not create new track from these hits " << endmsg ;
            return 0;            
          }
          
        }// for inner
      }// for outer    
    }// for middle
  }// if all vectors are not empty
  //    float dZ = FastTripletCheck(innerHit, middleHit, outerHit);
  
  //    if (fabs(dZ) > _minDistCutAttach)
  //      return trackAR;    
  
  // increase triplet count
  ++_ntriplets;

  // get the hit coordinates and errors
  double xh[3];
  double yh[3];
  float  zh[3];
  double wrh[3];
  float  wzh[3];
  float  rh[3];
  float  ph[3];
  
  float par[5];
  float epar[15];

  // first hit
  xh[0] = outerHit->getTrackerHit().getPosition()[0];
  yh[0] = outerHit->getTrackerHit().getPosition()[1];
  zh[0] = float(outerHit->getTrackerHit().getPosition()[2]);
  wrh[0] = double(1.0/(outerHit->getResolutionRPhi()*outerHit->getResolutionRPhi()));
  wzh[0] = 1.0/(outerHit->getResolutionZ()*outerHit->getResolutionZ());
  // second hit
  xh[1] = middleHit->getTrackerHit().getPosition()[0];
  yh[1] = middleHit->getTrackerHit().getPosition()[1];
  zh[1] = float(middleHit->getTrackerHit().getPosition()[2]);
  wrh[1] = double(1.0/(middleHit->getResolutionRPhi()*middleHit->getResolutionRPhi()));
  wzh[1] = 1.0/(middleHit->getResolutionZ()*middleHit->getResolutionZ());
  // third hit
  xh[2] = innerHit->getTrackerHit().getPosition()[0];
  yh[2] = innerHit->getTrackerHit().getPosition()[1];
  zh[2] = float(innerHit->getTrackerHit().getPosition()[2]);
  wrh[2] = double(1.0/(innerHit->getResolutionRPhi()*innerHit->getResolutionRPhi()));
  wzh[2] = 1.0/(innerHit->getResolutionZ()*innerHit->getResolutionZ());
  // calculate r and phi for all hits
  for (int ih=0; ih<3; ih++) {
    rh[ih] = float(sqrt(xh[ih]*xh[ih]+yh[ih]*yh[ih]));
    ph[ih] = atan2(yh[ih],xh[ih]);
    if (ph[ih] < 0.) 
      ph[ih] = TWOPI + ph[ih]; 
  }
  
  int NPT = 3;
  int iopt = 2;
  float chi2RPhi;
  float chi2Z;
  
  debug() << " TestTriplet: Use fastHelixFit " << endmsg ;  
  
  _fastfitter->fastHelixFit(NPT, xh, yh, rh, ph, wrh, zh, wzh,iopt, par, epar, chi2RPhi, chi2Z);
  par[3] = par[3]*par[0]/fabs(par[0]);

  // get helix parameters
  float omega = par[0];
  float tanlambda = par[1];
  float phi0 = par[2];
  float d0 = par[3];
  float z0 = par[4];

  // chi2 is weighted here by a factor for both rphi and z
  float Chi2 = chi2RPhi*_chi2WRPhiTriplet+chi2Z*_chi2WZTriplet;
  int ndf = 2*NPT-5;

  
  // check the truth information for the triplet
  
  // define these outside of the ifdef so that we don't need to keep repeating it.
  //std::vector<ConstTrackerHit> hit_list;
  //std::vector<MCParticle*> mcps_imo;
  //std::vector<MCParticle*> mcp_s;
  int triplet_code = 0;
  /*  
#ifdef MARLINTRK_DIAGNOSTICS_ON

  int nmcps   = 0;
  int nbadHits = 0;
  
  int layer  = 9 ;
  int size   = 3 ;
  int marker = 1 ;
  int ml     = 0 ;
  //  float helix_max_r = 0;
  float helix_max_z = 0;
  int color = 0;

  // use the MCTruth4HitExt to get the MCPs
  
  hit_list.push_back(innerHit->getTrackerHit());
  hit_list.push_back(middleHit->getTrackerHit());
  hit_list.push_back(outerHit->getTrackerHit());

  EVENT::MCParticle* mcp_i = 0;
  EVENT::MCParticle* mcp_m = 0;
  EVENT::MCParticle* mcp_o = 0;
  
  for (unsigned ihit = 0; ihit < hit_list.size(); ++ihit) {

    EVENT::ConstTrackerHit trkhit = hit_list[ihit];
    std::vector<MCParticle*> mcps;

    MarlinTrk::getMCParticlesForTrackerHit(trkhit, mcps);
    
    if (mcps.size() == 1) {
      mcps_imo.push_back(mcps[0]);
      ++nmcps;
    } else {
      mcps_imo.push_back(0);
      ++nbadHits;
    }
    
  }
    
  mcp_i = mcps_imo[0];
  mcp_m = mcps_imo[1];
  mcp_o = mcps_imo[2];
  
  debug()
  << "\n mcp_i = " << mcp_i
  << "\n mcp_m = " << mcp_m
  << "\n mcp_o = " << mcp_o
  << endmsg;
  
  if( mcp_i ) {
    mcp_s.push_back(mcp_i) ;
  }
    
  if( mcp_m && mcp_m != mcp_i ) {
    mcp_s.push_back(mcp_m);
  }
  
  if( mcp_o && mcp_o != mcp_m && mcp_o != mcp_i ){
    mcp_s.push_back(mcp_o);
  }

  nmcps = mcp_s.size();

  if (_UseEventDisplay) {
    // display this triplet and the MCPs from which it is formed
    
    MarlinCED::newEvent(this , _detector_model_for_drawing ) ;
    
    //    CEDPickingHandler &pHandler=CEDPickingHandler::getInstance();
    //
    //    pHandler.update(_current_event);
    
    for (unsigned imcp = 0; imcp < mcp_s.size(); ++imcp) {
      
      MCParticle* mcp = mcp_s[imcp];
      
      helix_max_z = fabsf(mcp->getEndpoint()[2]);
      
      
      info() << "Draw MCParticle : " << *mcp <<endmsg;
      
      MarlinCED::add_layer_description("MCParticle_For_Fit", layer);
      
      MarlinCED::drawHelix( _bField , mcp->getCharge(), mcp->getVertex()[0], mcp->getVertex()[1], mcp->getVertex()[2],
                           mcp->getMomentum()[0], mcp->getMomentum()[1], mcp->getMomentum()[2], layer , size , 0x7af774  ,
                           0.0,  _helix_max_r ,
                           helix_max_z, mcp.id() ) ;
      
    }
    
    const std::string  colName = "Hits_For_Fit";
    
    
    size   = 10 ;
    layer  = 11 ;
    //    ml = marker | ( layer << CED_LAYER_SHIFT ) ;
    
    //ced_describe_layer( colName.c_str() ,layer);
    MarlinCED::add_layer_description(colName, layer);
    
    
    color =  0xFFFFFF;
    
    for(   std::vector<ConstTrackerHit >::const_iterator it = hit_list.begin();  it != hit_list.end() ; it++ ) {
      
      ConstTrackerHit trkhit = *it;
      
      ced_hit_ID(trkhit.getPosition()[0],
                 trkhit.getPosition()[1],
                 trkhit.getPosition()[2],
                 marker, layer, size , color, trkhit.id() ) ;
      
    } // hits
    
  }
  
  if (_createDiagnosticsHistograms) {
    
    // if no bad hits are present triplet_code = nmcps;
    triplet_code = nmcps + nbadHits * 3  ;
        
    _histos->fill1D(DiagnosticsHistograms::htriplets, triplet_code);
    
    double pt =  (2.99792458E-4*_bField) / omega ; // for r in mm, p in GeV and Bz in Tesla
    
    if (triplet_code == 1) {
      ++_ntriplets_good;
      _histos->fill2D(DiagnosticsHistograms::htripletChi2vPt_good, pt, Chi2 );
      _histos->fill1D(DiagnosticsHistograms::htriplets_chi2_good, Chi2 );
      _histos->fill1D(DiagnosticsHistograms::htriplets_pt_good, pt );
    } else {

      _histos->fill2D(DiagnosticsHistograms::htripletChi2vPt_bad, pt, Chi2);
      _histos->fill1D(DiagnosticsHistograms::htriplets_chi2_bad, Chi2 );
      _histos->fill1D(DiagnosticsHistograms::htriplets_pt_bad, pt );

      if(triplet_code == 2) {
        ++_ntriplets_2MCP;
      } else if (triplet_code == 3) {
        ++_ntriplets_3MCP;
      } else if (triplet_code == 4) {
        ++_ntriplets_1MCP_Bad;
      } else {
        ++_ntriplets_bad;
      }
    }
    
  }
  
#endif
  */  
  
  // Check if track satisfies all conditions
  
  
  //   std::cout << "Chi2/ndf = " << Chi2/float(ndf) << " , cut = " << _chi2FitCut << endmsg;
  //   std::cout << "d0 = " << d0 << " , cut = " << _cutOnD0  << endmsg;
  //   std::cout << "z0 = " << z0 << " , cut = " << _cutOnZ0  << endmsg;
  //   std::cout << "omega = " << omega << " , cut = " << _cutOnOmega << endmsg;
  
  //  if ( Chi2/float(ndf) > _chi2FitCut || fabs(d0) > _cutOnD0 || fabs(z0) > _cutOnZ0 || fabs(omega)>_cutOnOmega)
  // return a null pointer
  //    return 0;
  
  bool failed = false;

  int quality_code = triplet_code * 10 ;

  if ( Chi2/float(ndf) > _chi2FitCut ) {
    debug() << "Chi2/ndf = " << Chi2/float(ndf) << " , cut = " << _chi2FitCut << endmsg;
    failed = true;
    quality_code += 1;
  } else if (fabs(d0) > _cutOnD0 ) {
    debug() << "d0 = " << d0 << " , cut = " << _cutOnD0  << endmsg;
    failed = true;
    quality_code += 2;
  } else if (fabs(z0) > _cutOnZ0 ) {
    debug() << "z0 = " << z0 << " , cut = " << _cutOnZ0  << endmsg;
    failed = true;
    quality_code += 3;
  } else if ( fabs(omega)>_cutOnOmega)  {
    debug() << "omega = " << omega << " , cut = " << _cutOnOmega << endmsg;
    failed = true;
    quality_code += 4;
  } else {
    debug() << "Success !!!!!!!" << endmsg;
  }
  /*
  if (_createDiagnosticsHistograms) _histos->fill1D(DiagnosticsHistograms::htriplets, quality_code);

    if (_UseEventDisplay) {
    drawEvent();
  }
  */
  
  if( failed ) {
    // return a null pointer
    return 0;
  }
  
  
  helix.Initialize_Canonical(phi0,d0,z0,omega,tanlambda,_bField);
  
  TrackExtended * trackAR = new TrackExtended();
  trackAR->addTrackerHitExtended(outerHit);
  trackAR->addTrackerHitExtended(middleHit);
  trackAR->addTrackerHitExtended(innerHit);
  outerHit->addTrackExtended(trackAR);
  middleHit->addTrackExtended(trackAR);
  innerHit->addTrackExtended(trackAR);    
  trackAR->setD0(d0);
  trackAR->setZ0(z0);
  trackAR->setPhi(phi0);
  trackAR->setTanLambda(tanlambda);
  trackAR->setOmega(omega);
  trackAR->setChi2( Chi2 );
  trackAR->setNDF( ndf );
  trackAR->setCovMatrix(epar);
  
  
  return trackAR;
  
}

int SiliconTrackingAlg::BuildTrack(TrackerHitExtended * outerHit, 
                                          TrackerHitExtended * middleHit,
                                          TrackerHitExtended * innerHit,
                                          HelixClass & helix,
                                          int innerLayer,
                                          int iPhiLow, int iPhiUp,
                                          int iThetaLow, int iThetaUp, 
                                          TrackExtended * trackAR) {
  /**
   Method for building up track in the VXD. Method starts from the found triplet and performs
   sequential attachment of hits in other layers, which have hits within the search window.
   Only searches inwards.
   Given that we know we are now jumping over layers due to the doublet nature of the VXD, we 
   could optimise this to look for the hits in interleaving layers as well. 
   Currently a fast fit is being done for each additional hit, it could be more efficient to try and use kaltest?
   
   */
  
  debug() << " BuildTrack starting " << endmsg;
  
  for (int layer = innerLayer-1; layer>=0; layer--) { // loop over remaining layers
    float distMin = 1.0e+20;
    TrackerHitExtended * assignedhit = NULL;
    // loop over phi in the Inner region
    for (int ipInner=iPhiLow; ipInner<iPhiUp+1;ipInner++) { 
      // loop over theta in the Inner region 
      for (int itInner=iThetaLow; itInner<iThetaUp+1;itInner++) { 
        int iPhiInner = ipInner;
        
        // catch wrap-around
        if (ipInner < 0) iPhiInner = _nDivisionsInPhi-1;
        if (ipInner >= _nDivisionsInPhi) iPhiInner = ipInner - _nDivisionsInPhi;
        
        // get the index of the theta-phi bin to search
        int iCode = layer + _nLayers*iPhiInner +  _nLayers*_nDivisionsInPhi*itInner;
        
        // get the hits from this bin
        TrackerHitExtendedVec& hitVecInner = _sectors[iCode];
        
        int nHitsInner = int(hitVecInner.size());
        
        // loop over hits in the Inner sector
        for (int iInner=0;iInner<nHitsInner;iInner++) { 
          
          TrackerHitExtended * currentHit = hitVecInner[iInner];
          
          // get the position of the hit to test
          float pos[3];
          float distance[3];
          
          for (int i=0; i<3; ++i) {
            pos[i] = float(currentHit->getTrackerHit().getPosition()[i]);
          }
          
          // get the distance of closest approach and distance s traversed to the POCA 
          float time = helix.getDistanceToPoint(pos,distance);    
          
          // sanity check on s 
          if (time < 1.0e+10) {
            
            // check if this is the closest hit yet
            if (distance[2] < distMin) { // distance[2] = sqrt( d0*d0 + z0*z0 ) 
              
              // if yes store hit and distance 
              distMin = distance[2];             
              assignedhit = currentHit;
            }
          }
        } // endloop over hits in the Inner sector
      } // endloop over theta in the Inner region 
    } // endloop over phi in the Inner region
    // check if closest hit fulfills the min distance cut
    if (distMin < _minDistCutAttach) {
      // if yes try to include it in the fit 
      TrackerHitExtendedVec& hvec = trackAR->getTrackerHitExtendedVec();
      int  nHits = int(hvec.size());
      double * xh = new double[nHits+1];
      double * yh = new double[nHits+1];
      float * zh = new float[nHits+1];
      double * wrh = new double[nHits+1];
      float * wzh = new float[nHits+1];
      float * rh = new float[nHits+1];
      float * ph = new float[nHits+1];
      float par[5];
      float epar[15];
      
      for (int ih=0;ih<nHits;++ih) {
	edm4hep::ConstTrackerHit trkHit = hvec[ih]->getTrackerHit();
        xh[ih] = trkHit.getPosition()[0];
        yh[ih] = trkHit.getPosition()[1];
        zh[ih] = float(trkHit.getPosition()[2]);
        wrh[ih] = double(1.0/(hvec[ih]->getResolutionRPhi()*hvec[ih]->getResolutionRPhi()));
        wzh[ih] = 1.0/(hvec[ih]->getResolutionZ()*hvec[ih]->getResolutionZ());
        rh[ih] = float(sqrt(xh[ih]*xh[ih]+yh[ih]*yh[ih]));
        ph[ih] = float(atan2(yh[ih],xh[ih]));
        if (ph[ih] < 0.) 
          ph[ih] = TWOPI + ph[ih]; 
      }      
      edm4hep::ConstTrackerHit assignedTrkHit = assignedhit->getTrackerHit();
      xh[nHits] = assignedTrkHit.getPosition()[0];
      yh[nHits] = assignedTrkHit.getPosition()[1];
      zh[nHits] = float(assignedTrkHit.getPosition()[2]);
      rh[nHits] = float(sqrt(xh[nHits]*xh[nHits]+yh[nHits]*yh[nHits]));
      ph[nHits] = float(atan2(yh[nHits],xh[nHits]));
      if (ph[nHits] < 0.) 
        ph[nHits] = TWOPI + ph[nHits]; 
      wrh[nHits] = double(1.0/(assignedhit->getResolutionRPhi()*assignedhit->getResolutionRPhi()));
      wzh[nHits] = 1.0/(assignedhit->getResolutionZ()*assignedhit->getResolutionZ());
      
      int NPT = nHits + 1;
      int iopt = 2;
      float chi2RPhi;
      float chi2Z;
      
      //debug() << "######## number of hits to fit with _fastfitter = " << NPT << endmsg; 
      
      _fastfitter->fastHelixFit(NPT, xh, yh, rh, ph, wrh, zh, wzh,iopt, par, epar, chi2RPhi, chi2Z);
      par[3] = par[3]*par[0]/fabs(par[0]);
      
      
      delete[] xh;
      delete[] yh;
      delete[] zh;
      delete[] wrh;
      delete[] wzh;
      delete[] rh;
      delete[] ph;
      
      bool validCombination = 0;
      float Chi2 = FLT_MAX;
      
      if ((nHits+1) == 4) {
        Chi2 = chi2RPhi*_chi2WRPhiQuartet+chi2Z*_chi2WZQuartet;
      }         
      if ((nHits+1) >= 5) {
        Chi2 = chi2RPhi*_chi2WRPhiSeptet+chi2Z*_chi2WZSeptet;
      }
      int ndf = 2*NPT-5;
      
      // check if this is valid combination based on the chi2/ndf
      validCombination = Chi2/float(ndf) < _chi2FitCut;
      
      if ( validCombination ) {
        // assign hit to track and track to hit, update the track parameters
        trackAR->addTrackerHitExtended(assignedhit);
        assignedhit->addTrackExtended(trackAR);
        float omega = par[0];
        float tanlambda = par[1];
        float phi0 = par[2];
        float d0 = par[3];
        float z0 = par[4];
        helix.Initialize_Canonical(phi0,d0,z0,omega,tanlambda,_bField);
        trackAR->setD0(d0);
        trackAR->setZ0(z0);
        trackAR->setPhi(phi0);
        trackAR->setTanLambda(tanlambda);
        trackAR->setOmega(omega);
        trackAR->setChi2( Chi2 );
        trackAR->setCovMatrix(epar);
        trackAR->setNDF( ndf );
      }
      
    }
  } // endloop over remaining layers
  TrackerHitExtendedVec& hvec = trackAR->getTrackerHitExtendedVec();  
  int nTotalHits = int(hvec.size());
  debug() << "######## number of hits to return = " << nTotalHits << endmsg; 
  return nTotalHits;
}


void SiliconTrackingAlg::Sorting(TrackExtendedVec & trackVec) {
  /**
   Sorting of Track Vector in ascending order of chi2/ndf
   */
  
  std::sort(trackVec.begin(), trackVec.end(), compare_TrackExtended() );
  
  // also clean up? what does this do here?
  for (size_t i=0, sizeOfVector=trackVec.size(); i<sizeOfVector; ++i) {
    
    TrackerHitExtendedVec& hitVec = trackVec[i]->getTrackerHitExtendedVec();
    int nHits = int(hitVec.size());
    
    for (int ih=0;ih<nHits;ih++) {
      hitVec[ih]->clearTrackVec();
    }
  }
  
}

void SiliconTrackingAlg::CreateTrack(TrackExtended * trackAR ) {
  
  /**
   Method which creates Track out of TrackExtended objects. Checks for possible
   track splitting (separate track segments in VXD and FTD).
   */
  
  
  TrackerHitExtendedVec& hitVec = trackAR->getTrackerHitExtendedVec();
  int nHits = int(hitVec.size());
  
  for (int i=0; i<nHits; ++i) {
    TrackExtendedVec& trackVec = hitVec[i]->getTrackExtendedVec();
    if (trackVec.size() != 0) 
      return ;
  }
  
  // First check if the current track is piece of the split one
  // look for matching track segment
  
  int found = 0;
  
  int nTrk = int(_trackImplVec.size());
  
  for (int itrk=0; itrk<nTrk; ++itrk) {
    TrackExtended * trackOld = _trackImplVec[itrk];
    TrackerHitExtendedVec& hitVecOld = trackOld->getTrackerHitExtendedVec();
    
    float phiNew = trackAR->getPhi();
    float phiOld = trackOld->getPhi();
    float thetaNew = M_PI_2 - atan(trackAR->getTanLambda());
    float thetaOld = M_PI_2 - atan(trackOld->getTanLambda());
    
    float angle = (cos(phiNew)*cos(phiOld)+sin(phiNew)*sin(phiOld))*sin(thetaNew)*sin(thetaOld)+cos(thetaNew)*cos(thetaOld);
    angle = acos(angle);
    
    if (angle < _angleCutForMerging) {
      int nHitsOld = int(hitVecOld.size());
      int nTotHits = nHits + nHitsOld;
      double * xh = new double[nTotHits];
      double * yh = new double[nTotHits];
      float * zh = new float[nTotHits];
      double * wrh = new double[nTotHits];
      float * wzh = new float[nTotHits];
      float * rh = new float[nTotHits];
      float * ph = new float[nTotHits];
      float par[5];
      float epar[15];
      float refPoint[3] = {0.,0.,0.};
      for (int ih=0;ih<nHits;++ih) {
	edm4hep::ConstTrackerHit trkHit = hitVec[ih]->getTrackerHit();
        float rR = hitVec[ih]->getResolutionRPhi();
        float rZ = hitVec[ih]->getResolutionZ();
        if (int(hitVec[ih]->getTrackExtendedVec().size()) != 0)
          debug() << "WARNING : HIT POINTS TO TRACK " << endmsg;
        xh[ih] = trkHit.getPosition()[0];
        yh[ih] = trkHit.getPosition()[1];
        zh[ih] = float(trkHit.getPosition()[2]);
        wrh[ih] = double(1.0/(rR*rR));
        wzh[ih] = 1.0/(rZ*rZ);
        rh[ih] = float(sqrt(xh[ih]*xh[ih]+yh[ih]*yh[ih]));
        ph[ih] = float(atan2(yh[ih],xh[ih]));
      }      
      for (int ih=0;ih<nHitsOld;++ih) {
	edm4hep::ConstTrackerHit trkHit = hitVecOld[ih]->getTrackerHit();
        xh[ih+nHits] = trkHit.getPosition()[0];
        yh[ih+nHits] = trkHit.getPosition()[1];
        zh[ih+nHits] = float(trkHit.getPosition()[2]);
        float rR = hitVecOld[ih]->getResolutionRPhi();
        float rZ = hitVecOld[ih]->getResolutionZ();     
        wrh[ih+nHits] = double(1.0/(rR*rR));
        wzh[ih+nHits] = 1.0/(rZ*rZ);
        rh[ih+nHits] = float(sqrt(xh[ih+nHits]*xh[ih+nHits]+yh[ih+nHits]*yh[ih+nHits]));
        ph[ih+nHits] = float(atan2(yh[ih+nHits],xh[ih+nHits]));
        
      }
      int NPT = nTotHits;
      int iopt = 2;
      float chi2RPhi;
      float chi2Z;
      int ndf = 2*NPT - 5;
      
      _fastfitter->fastHelixFit(NPT, xh, yh, rh, ph, wrh, zh, wzh,iopt, par, epar, chi2RPhi, chi2Z);
      par[3] = par[3]*par[0]/fabs(par[0]);
      
      float omega = par[0];
      float tanlambda = par[1];
      float phi0 = par[2];
      float d0 = par[3];
      float z0 = par[4];
      
      float eparmin[15];
      for (int iparam=0;iparam<15;++iparam)
        eparmin[iparam] = epar[iparam];      
      
      float refPointMin[3];
      for (int ipp=0;ipp<3;++ipp)
        refPointMin[ipp] = refPoint[ipp];
      
      float chi2Min = chi2RPhi*_chi2WRPhiSeptet+chi2Z*_chi2WZSeptet;
      chi2Min = chi2Min/float(ndf);
      
      float chi2MinRPhi = chi2RPhi;
      float chi2MinZ = chi2Z;
      
      int iBad = -1;
      if (chi2Min < _chi2FitCut) {
        found = 1;
      }
      else { // SJA:FIXME: UH What is going on here? setting weights to 0 and refitting?
        float * wzhOld = new float[nTotHits];
        double * wrhOld = new double[nTotHits];
        for (int i=0;i<nTotHits;++i) {
          wzhOld[i] = wzh[i];
          wrhOld[i] = wrh[i];
        }
        for (int i=0; i<nTotHits; ++i) {
          for (int j=0;j<nTotHits;++j) {
            if (i == j) {
              wrh[j] = 0.0;
              wzh[j] = 0.0;
            } 
            else {
              wrh[j] = wrhOld[j];
              wzh[j] = wzhOld[j];
            }
          }
          
          _fastfitter->fastHelixFit(NPT, xh, yh, rh, ph, wrh, zh, wzh,iopt, par, epar, chi2RPhi, chi2Z);
          par[3] = par[3]*par[0]/fabs(par[0]);
          
          float chi2Cur = chi2RPhi*_chi2WRPhiSeptet+chi2Z*_chi2WZSeptet;
          chi2Cur = chi2Cur/float(ndf);
          
          if (chi2Cur < chi2Min) {
            chi2Min = chi2Cur;
            chi2MinRPhi = chi2RPhi;
            chi2MinZ = chi2Z;
            omega = par[0];
            tanlambda = par[1];
            phi0 = par[2];
            d0 = par[3];
            z0 = par[4];
            for (int iparam=0;iparam<15;++iparam)
              eparmin[iparam] = epar[iparam];
            for (int ipp=0;ipp<3;++ipp)
              refPointMin[ipp] = refPoint[ipp];
            iBad = i;
          }
        }
        if (chi2Min < _chi2FitCut) {
          found = 1;
        }
        delete[] wzhOld;
        delete[] wrhOld;
      }
      
      // Split track is found.
      // Attach hits belonging to the current track segment to  
      // the track already created
      if (found == 1) {
        trackOld->ClearTrackerHitExtendedVec();
        for (int i=0;i<nHits;++i) {
          TrackerHitExtended * trkHit = hitVec[i];
          trkHit->clearTrackVec();
          if (i == iBad) {          
          }
          else {
            trackOld->addTrackerHitExtended(trkHit);
            trkHit->addTrackExtended( trackOld );
          }
        }  
        for (int i=0;i<nHitsOld;++i) {
          int icur = i+nHits;
          TrackerHitExtended * trkHit = hitVecOld[i];
          trkHit->clearTrackVec();
          if (icur == iBad) {
          }
          else {
            trackOld->addTrackerHitExtended(trkHit);
            trkHit->addTrackExtended( trackOld );
          }
        }
        trackOld->setOmega(omega);
        trackOld->setTanLambda(tanlambda);
        trackOld->setPhi(phi0);
        trackOld->setD0(d0);
        trackOld->setZ0(z0);
        
        //debug() << "Split track found " << d0 << " " << z0 << endmsg;
        
        // killeb:  In the original SiliconTrackingAlg this was in the NOT simple helix branch.
        // The rest of the code uses the simple helix branch, where ndf_D is never set.
        // In fact it has never been initialised or used anywhere. I think this line should not be executed.
        // ndf = ndf_D;
        
        trackOld->setChi2(chi2Min*float(ndf));  
        trackOld->setNDF(ndf);
        trackOld->setCovMatrix(eparmin);
        //      trackOld->setReferencePoint(refPointMin);
      }
      
      delete[] xh;
      delete[] yh;
      delete[] zh;
      delete[] wrh;
      delete[] wzh;
      delete[] rh;
      delete[] ph;
      
    }
    if (found == 1)
      break;
  }
  
  // Candidate is a unique track
  // No other segments are found
  if (found == 0 ) {
    _trackImplVec.push_back(trackAR);
    for (int i=0;i<nHits;++i) {
      TrackerHitExtended * hit = hitVec[i];
      hit->addTrackExtended( trackAR );
    }
  }
  
  
}

void SiliconTrackingAlg::AttachRemainingVTXHitsFast() {
  
  std::vector<TrackerHitExtendedVec> nonAttachedHits;
  nonAttachedHits.resize(_nDivisionsInPhi*_nDivisionsInTheta);
  std::vector<TrackExtendedVec> trackVector;
  trackVector.resize(_nDivisionsInPhi*_nDivisionsInTheta);
  int nTracks = int(_trackImplVec.size());
  
  for (int iTrk=0;iTrk<nTracks;++iTrk) {
    TrackExtended * track = _trackImplVec[iTrk];
    double Phi = double(track->getPhi());
    if (Phi < 0)
      Phi = Phi + TWOPI;
    float tanlambda = track->getTanLambda();
    double cosTheta = double(tanlambda/sqrt(1+tanlambda*tanlambda));
    int iPhi = int(Phi/_dPhi);
    int iTheta = int ((cosTheta + double(1.0))/_dTheta);
    int iCode = iPhi + _nDivisionsInPhi*iTheta; 
    trackVector[iCode].push_back( track );
  }
  
  for (int il=0;il<_nLayers;++il) {
    for (int ip=0;ip<_nDivisionsInPhi;++ip) {
      for (int it=0;it<_nDivisionsInTheta; ++it) {
        int iCode = il + _nLayers*ip + _nLayers*_nDivisionsInPhi*it;      
        TrackerHitExtendedVec& hitVec = _sectors[iCode];
        int nH = int(hitVec.size());
        for (int iH=0; iH<nH; ++iH) {
          TrackerHitExtended * hitExt = hitVec[iH];
          TrackExtendedVec& trackVec = hitExt->getTrackExtendedVec();
          if (trackVec.size()==0) {
	    edm4hep::ConstTrackerHit hit = hitExt->getTrackerHit();
            double pos[3];
            double radius = 0;
            for (int i=0; i<3; ++i) {
              pos[i] = hit.getPosition()[i];
              radius += pos[i]*pos[i];
            }
            radius = sqrt(radius);
            double cosTheta = pos[2]/radius;
            double Phi = atan2(pos[1],pos[0]);
            if (Phi < 0.) Phi = Phi + TWOPI;
            int iPhi = int(Phi/_dPhi);
            int iTheta = int ((cosTheta + double(1.0))/_dTheta);
            iCode = iPhi + _nDivisionsInPhi*iTheta;      
            nonAttachedHits[iCode].push_back( hitExt );
          }
        }
      }
    }
  }
  
  for (int iT=0; iT<_nDivisionsInTheta; ++iT) {
    for (int iP=0; iP<_nDivisionsInPhi; ++iP) {
      int iCode = iP + _nDivisionsInPhi*iT; 
      int nHits = int(nonAttachedHits[iCode].size());
      int iT1 = iT - 1;
      int iT2 = iT + 1; 
      if (iT == 0) {
        iT1 = iT;
        iT2 = iT1 + 1;
      }
      if (iT == _nDivisionsInTheta - 1) {
        iT2 = iT;
        iT1 = iT2 - 1;
      }
      int iPHI[3];
      iPHI[0] = iP - 1;
      iPHI[1] = iP;
      iPHI[2] = iP + 1;
      if (iP == 0) 
        iPHI[0] = _nDivisionsInPhi - 1;
      if (iP == _nDivisionsInPhi - 1 )
        iPHI[2] = 0;
      
      for (int ihit = 0; ihit<nHits; ++ihit) {
        
        TrackerHitExtended * hit = nonAttachedHits[iCode][ihit];
        TrackExtended * trackToAttach = NULL;
        float minDist = 1.0e+6;
        
        for (int iTheta = iT1; iTheta <iT2+1; ++iTheta) {
          for (int indexP=0;indexP<3;++indexP) {
            int iPhi = iPHI[indexP];        
            int iCodeForTrack = iPhi + _nDivisionsInPhi*iTheta;
            int nTrk = int(trackVector[iCodeForTrack].size());
            for (int iTrk=0; iTrk<nTrk; ++iTrk) {         
              TrackExtended * trackAR = trackVector[iCodeForTrack][iTrk];
              bool consider = true;
              if (_checkForDelta) {
                TrackerHitExtendedVec& hitVector = trackAR->getTrackerHitExtendedVec();
                int NHITS = int(hitVector.size());
                for (int IHIT=0;IHIT<NHITS;++IHIT) {
                  
                  // Here we are trying to find if a hits are too close i.e. closer than _minDistToDelta
		  edm4hep::ConstTrackerHit trkhit1 = hit->getTrackerHit();
		  edm4hep::ConstTrackerHit trkhit2 = hitVector[IHIT]->getTrackerHit();                  
                  
                  if ( trkhit1.getCellID() == trkhit2.getCellID() ){ // i.e. they are in the same sensor
                    float distance = 0.;
                    for (int iC=0;iC<3;++iC) {
                      float posFirst = float(hit->getTrackerHit().getPosition()[iC]);
                      float posSecond = float(hitVector[IHIT]->getTrackerHit().getPosition()[iC]);
                      float deltaPos = posFirst - posSecond;
                      distance += deltaPos*deltaPos;
                    }
                    distance = sqrt(distance);
                    if (distance<_minDistToDelta) {
                      consider = false;
                      break;
                    }
                  }
                }
              }
              if (consider) {   
                float phi0 = trackAR->getPhi();
                float d0 = trackAR->getD0();
                float z0 = trackAR->getZ0();
                float omega = trackAR->getOmega();
                float tanlambda = trackAR->getTanLambda();
                HelixClass helix;
                helix.Initialize_Canonical(phi0,d0,z0,omega,tanlambda,_bField);
                int layer = getLayerID(hit->getTrackerHit());
                if (layer > _minimalLayerToAttach) {
                  float pos[3];
                  for (int i=0; i<3; ++i) 
                    pos[i] = hit->getTrackerHit().getPosition()[i];      
                  float distance[3];
                  float time = helix.getDistanceToPoint(pos,distance);
                  if (time < 1.0e+10) {
                    if (distance[2] < minDist) {
                      minDist = distance[2];
                      trackToAttach = trackAR;
                    }                      
                  }    
                }
              }
            }
          }
        }
        if (minDist < _minDistCutAttach && trackToAttach != NULL) {
          int iopt = 2;
          AttachHitToTrack(trackToAttach,hit,iopt);
        }      
      }
    }
  }
}

void SiliconTrackingAlg::AttachRemainingVTXHitsSlow() {
  
  TrackerHitExtendedVec nonAttachedHits;
  nonAttachedHits.clear();
  
  for (int il=0;il<_nLayers;++il) {
    for (int ip=0;ip<_nDivisionsInPhi;++ip) {
      for (int it=0;it<_nDivisionsInTheta; ++it) {
        int iCode = il + _nLayers*ip + _nLayers*_nDivisionsInPhi*it;      
        TrackerHitExtendedVec& hitVec = _sectors[iCode];
        int nH = int(hitVec.size());
        for (int iH=0; iH<nH; ++iH) {
          TrackerHitExtended * hit = hitVec[iH];
          TrackExtendedVec& trackVec = hit->getTrackExtendedVec();
          // if (trackVec.size()==0)
          // nonAttachedHits.push_back( hit );
          //-- allow hits that are only used in triplets to be re-attached 
          unsigned int maxTrackSize = 0;
          for(unsigned int itrack = 0;itrack < trackVec.size();itrack++){
            TrackerHitExtendedVec hitVec_tmp= trackVec[itrack]->getTrackerHitExtendedVec();
            unsigned int isize = hitVec_tmp.size();
            if(isize>maxTrackSize)maxTrackSize = isize;
          }     
          if (maxTrackSize<=3) { 
            debug() << " Add non attached hit to list: id = " << hit->getTrackerHit().id() << endmsg;
            nonAttachedHits.push_back( hit );
          } 
          
          
        }
      }
    }
  }
  
  int nNotAttached = int(nonAttachedHits.size());
  
  int nTrk = int(_trackImplVec.size()); 
  for (int iHit=0; iHit<nNotAttached; ++iHit) {
    TrackerHitExtended * hit = nonAttachedHits[iHit];
    debug() << " Try hit: id = " << hit->getTrackerHit().id() << endmsg;
    int layer = getLayerID( hit->getTrackerHit() );
    if (layer > _minimalLayerToAttach) {
      float pos[3];
      for (int i=0; i<3; ++i) 
        pos[i] = hit->getTrackerHit().getPosition()[i];      
      float minDist = 1e+10;
      TrackExtended * trackToAttach = NULL;
      for (int iTrk=0; iTrk<nTrk; ++iTrk) {
        TrackExtended * trackAR = _trackImplVec[iTrk];
        bool consider = true;
        if (_checkForDelta) {
          TrackerHitExtendedVec& hitVector = trackAR->getTrackerHitExtendedVec();
          int NHITS = int(hitVector.size());
          for (int IHIT=0;IHIT<NHITS;++IHIT) {
            
            // Here we are trying to find if a hits are too close i.e. closer than _minDistToDelta
	    edm4hep::ConstTrackerHit trkhit1 = hit->getTrackerHit();
	    edm4hep::ConstTrackerHit trkhit2 = hitVector[IHIT]->getTrackerHit();                  
	    
            if ( trkhit1.getCellID() == trkhit2.getCellID() ){ // i.e. they are in the same sensor
              
              float distance = 0.;
              for (int iC=0;iC<3;++iC) {
                float posFirst = float(hit->getTrackerHit().getPosition()[iC]);
                float posSecond = float(hitVector[IHIT]->getTrackerHit().getPosition()[iC]);
                float deltaPos = posFirst - posSecond;
                distance += deltaPos*deltaPos;
              }
              distance = sqrt(distance);
              if (distance<_minDistToDelta) {
                consider = false;
                debug() << " hit: id = " << hit->getTrackerHit().id() << " condsidered delta together with hit " << trkhit2.id() << endmsg;
                break;
              }
            }       
          }
        }
        if (consider) {
          HelixClass helix;
          float phi0 = trackAR->getPhi();
          float d0 = trackAR->getD0();
          float z0 = trackAR->getZ0();
          float omega = trackAR->getOmega();
          float tanlambda = trackAR->getTanLambda();
          helix.Initialize_Canonical(phi0,d0,z0,omega,tanlambda,_bField);
          float distance[3];
          float time = helix.getDistanceToPoint(pos,distance);
          if (time < 1.0e+10) {
            if (distance[2] < minDist) {
              minDist = distance[2];
              trackToAttach = trackAR;
            }
          }
        }
      }
      if (minDist < _minDistCutAttach && trackToAttach != NULL) {
        int iopt = 2;
        debug() << " Hit: id = " << hit->getTrackerHit().id() << " : try attachement"<< endmsg;
        AttachHitToTrack(trackToAttach,hit,iopt);
      } else {
        debug() << " Hit: id = " << hit->getTrackerHit().id() << " rejected due to distance cut of " <<_minDistCutAttach<< " min distance = "  << minDist << endmsg;
      }      
    }
  }  
}

void SiliconTrackingAlg::AttachRemainingFTDHitsSlow() {
  TrackerHitExtendedVec nonAttachedHits;
  nonAttachedHits.clear();
  
  for (int iS=0;iS<2;++iS) {
    for (unsigned int layer=0;layer<_nlayersFTD;++layer) {
      for (int ip=0;ip<_nDivisionsInPhiFTD;++ip) {
        int iCode = iS + 2*layer + 2*_nlayersFTD*ip;      
        TrackerHitExtendedVec& hitVec = _sectorsFTD[iCode];
        int nH = int(hitVec.size());
        for (int iH=0; iH<nH; ++iH) {
          TrackerHitExtended * hit = hitVec[iH];
          TrackExtendedVec& trackVec = hit->getTrackExtendedVec();
          if (trackVec.size()==0)
            nonAttachedHits.push_back( hit );
        }
      }
    }
  }
  
  int nNotAttached = int(nonAttachedHits.size());
  
  int nTrk = int(_trackImplVec.size()); 
  for (int iHit=0; iHit<nNotAttached; ++iHit) {
    TrackerHitExtended * hit = nonAttachedHits[iHit];
    float pos[3];
    for (int i=0; i<3; ++i) 
      pos[i] = hit->getTrackerHit().getPosition()[i];      
    float minDist = 1e+10;
    TrackExtended * trackToAttach = NULL;
    for (int iTrk=0; iTrk<nTrk; ++iTrk) {
      TrackExtended * trackAR = _trackImplVec[iTrk];
      bool consider = true;
      TrackerHitExtendedVec& hitVector = trackAR->getTrackerHitExtendedVec();
      int NHITS = int(hitVector.size());
      
      for (int IHIT=0;IHIT<NHITS;++IHIT) {
        
        // SJA:FIXME: check to see if allowing no hits in the same sensor vs no hits in the same layer works 
        //        if (hit->getTrackerHit()->getType() == hitVector[IHIT]->getTrackerHit()->getType()) {
        if (hit->getTrackerHit().getCellID() == hitVector[IHIT]->getTrackerHit().getCellID()) {
          
          consider = false;
          break;
        }
      }
      
      
      if (consider) {
        HelixClass helix;
        float phi0 = trackAR->getPhi();
        float d0 = trackAR->getD0();
        float z0 = trackAR->getZ0();
        float omega = trackAR->getOmega();
        float tanlambda = trackAR->getTanLambda();
        if (tanlambda*float(getSideID(hit->getTrackerHit())) > 0) {
          helix.Initialize_Canonical(phi0,d0,z0,omega,tanlambda,_bField);
          float distance[3];
          float time = helix.getDistanceToPoint(pos,distance);
          if (time < 1.0e+10) {
            if (distance[2] < minDist) {
              minDist = distance[2];
              trackToAttach = trackAR;
            }
          }
        }
      }
    }
    if (minDist < _minDistCutAttach && trackToAttach != NULL) {
      int iopt = 2;
      AttachHitToTrack(trackToAttach,hit,iopt);
    }      
  }  
}


void SiliconTrackingAlg::AttachRemainingFTDHitsFast() {
  int nTrk = _trackImplVec.size();
  
  for (int iTrk=0; iTrk<nTrk; ++iTrk) {
    TrackExtended * trackAR = _trackImplVec[iTrk];
    HelixClass helix;
    float phi0 = trackAR->getPhi();
    float d0 = trackAR->getD0();
    float z0 = trackAR->getZ0();
    float omega = trackAR->getOmega();
    float tanlambda = trackAR->getTanLambda();
    helix.Initialize_Canonical(phi0,d0,z0,omega,tanlambda,_bField);
    int iSemiSphere = 0;
    if (tanlambda > 0) 
      iSemiSphere = 1;
    float ref[3];
    for (int i=0;i<3;++i) 
      ref[i] = helix.getReferencePoint()[i];
    // Start loop over FTD layers
    for (unsigned int layer=0;layer<_nlayersFTD;layer++) {
      float ZL = _zLayerFTD[layer];
      if (iSemiSphere == 0)
        ZL = - ZL;
      float point[3];
      helix.getPointInZ(ZL,ref,point);
      float Phi = atan2(point[1],point[0]);
      if (Phi < 0) 
        Phi = Phi + TWOPI;
      int iPhi = int(Phi/_dPhiFTD);
      float distMin = 1e+10;
      TrackerHitExtended * attachedHit = NULL;     
      for (int iP=iPhi-1;iP<=iPhi+1;++iP) {
        int iPP = iP;
        if (iP < 0) 
          iPP = iP + _nDivisionsInPhiFTD;
        if (iP >= _nDivisionsInPhiFTD)
          iPP = iP - _nDivisionsInPhiFTD;  
        int iCode = iSemiSphere + 2*layer + 2*_nlayersFTD*iPP;
        int nHits = int(_sectorsFTD[iCode].size());
        for (int iHit=0;iHit<nHits;++iHit) {
          TrackerHitExtended * hit = _sectorsFTD[iCode][iHit];
          bool consider = true;
          TrackerHitExtendedVec& hitVector = trackAR->getTrackerHitExtendedVec();
          int NHITS = int(hitVector.size());
          
          // SJA:FIXME: check to see if allowing no hits in the same sensor vs no hits in the same layer works 
          for (int IHIT=0;IHIT<NHITS;++IHIT) {
            //            if (hit->getTrackerHit()->getType() == hitVector[IHIT]->getTrackerHit()->getType()) {
            if (hit->getTrackerHit().getCellID() == hitVector[IHIT]->getTrackerHit().getCellID()) {
              consider = false;
              break;
            }
          }
          
          
          if (consider) {
            float pos[3];
            for (int i=0;i<3;++i) {
              pos[i] = hit->getTrackerHit().getPosition()[i];
            }
            float distance[3];
            float time = helix.getDistanceToPoint(pos,distance);
            if (time < 1.0e+10) {
              if (distance[2] < distMin) {
                distMin = distance[2];
                attachedHit = hit;
              }
            }
          }
        }
      }
      if (distMin < _minDistCutAttach && attachedHit != NULL) {
        int iopt = 2;
        AttachHitToTrack(trackAR,attachedHit, iopt);
      }
    }
  }
}

void SiliconTrackingAlg::TrackingInFTD() {

  int nComb = int(_CombinationsFTD.size()) / 3;

  for (int iComb=0;iComb<nComb;++iComb) {

    int nLS[3];
    nLS[0] = _CombinationsFTD[3*iComb];
    nLS[1] = _CombinationsFTD[3*iComb+1];
    nLS[2] = _CombinationsFTD[3*iComb+2];

    for (int iS=0;iS<2;++iS) { // loop over +z and -z
      
      //      std::cout << "Combinations : " << iS << " " << nLS[0] << " " << nLS[1] << " " << nLS[2] << endmsg;
      //      int iC = iS + 2*nLS[0];
      //      TrackerHitExtendedVec& hitVec = _sectorsFTD[iC];
      //      int nO = int(hitVec.size());
      //      iC = iS + 2*nLS[1];
      //      hitVec = _sectorsFTD[iC];
      //      int nM = int(hitVec.size());
      //      iC = iS + 2*nLS[2];
      //      hitVec = _sectorsFTD[iC];
      //      int nI = int(hitVec.size());
      //      std::cout << nO << " " << nM << " " << nI << endmsg;

      for (int ipOuter=0;ipOuter<_nDivisionsInPhiFTD;++ipOuter) { 

        int ipMiddleLow = ipOuter - 1;
        int ipMiddleUp  = ipOuter + 1;
        
        unsigned int iCodeOuter = iS + 2*nLS[0] + 2*_nlayersFTD*ipOuter;
        
        if( iCodeOuter >= _sectorsFTD.size()){          
          error() << "iCodeOuter index out of range: iCodeOuter =   " << iCodeOuter << " _sectorsFTD.size() = " << _sectorsFTD.size() << " exit(1) called from file " << __FILE__ << " line " << __LINE__<< endmsg;
          exit(1);
        }
        
        TrackerHitExtendedVec& hitVecOuter = _sectorsFTD[iCodeOuter];
        
        int nOuter = int(hitVecOuter.size());

        for (int iOuter=0;iOuter<nOuter;++iOuter) {
          
          TrackerHitExtended * hitOuter = hitVecOuter[iOuter];
        
          for (int ipMiddle=ipMiddleLow;ipMiddle<=ipMiddleUp;++ipMiddle) {
            //for(int ipMiddle=0;ipMiddle<_nDivisionsInPhiFTD;++ipMiddle) {
            int ipM = ipMiddle;
            if (ipM < 0) 
              ipM = ipMiddle + _nDivisionsInPhiFTD;
            if (ipM >= _nDivisionsInPhiFTD)
              ipM = ipMiddle - _nDivisionsInPhiFTD;
            int iCodeMiddle = iS + 2*nLS[1] + 2*_nlayersFTD*ipM;
          
            TrackerHitExtendedVec& hitVecMiddle = _sectorsFTD[iCodeMiddle];
            int ipInnerLow,ipInnerUp;       
            ipInnerLow = ipMiddle - 1;
            ipInnerUp =  ipMiddle + 1;
            
            int nMiddle = int(hitVecMiddle.size());
            
            for (int iMiddle=0;iMiddle<nMiddle;++iMiddle) {
              TrackerHitExtended * hitMiddle = hitVecMiddle[iMiddle];
              for (int ipInner=ipInnerLow;ipInner<=ipInnerUp;++ipInner) {
                //for (int ipInner=0;ipInner<_nDivisionsInPhiFTD;++ipInner) {
                int ipI = ipInner;
                if (ipI < 0)
                  ipI = ipInner + _nDivisionsInPhiFTD;
                if (ipI >= _nDivisionsInPhiFTD) 
                  ipI = ipInner - _nDivisionsInPhiFTD;
                int iCodeInner = iS + 2*nLS[2] + 2*_nlayersFTD*ipI;
                TrackerHitExtendedVec& hitVecInner = _sectorsFTD[iCodeInner];
            
                int nInner = int(hitVecInner.size());
                
                for (int iInner=0;iInner<nInner;++iInner) {
                
                  TrackerHitExtended * hitInner = hitVecInner[iInner];
                  HelixClass helix;
                  //                  std::cout << endmsg;
                  //                  std::cout << "Outer Hit Type " << hitOuter->getTrackerHit()->getType() << " z = " << hitOuter->getTrackerHit().getPosition()[2] 
                  //                  << "\nMiddle Hit Type "<< hitMiddle->getTrackerHit()->getType() << " z = " << hitMiddle->getTrackerHit().getPosition()[2]  
                  //                  << "\nInner Hit Type "<< hitInner->getTrackerHit()->getType() << " z = " << hitInner->getTrackerHit().getPosition()[2]  << endmsg;

                  debug() << " "
                  << std::setw(3) << ipOuter       << " "   << std::setw(3) << ipMiddle << " "      << std::setw(3) << ipInner << "       "
                  << std::setw(3) << iS << "     "
                  << std::setw(3) << nLS[0]     << " "   << std::setw(3) << nLS[1]   << " "      << std::setw(3) << nLS[2]  << "     "
                  << std::setw(3) << nOuter << " : " << std::setw(3) << nMiddle << " : " << std::setw(3) << nInner << "  :: "
                  << std::setw(3) << nOuter*nMiddle* nInner << endmsg;

                  
                  TrackExtended * trackAR = TestTriplet(hitOuter,hitMiddle,hitInner,helix);
                  if (trackAR != NULL) {
                    //                    std::cout << "FTD triplet found" << endmsg;
                    int nHits = BuildTrackFTD(trackAR,nLS,iS);
                    
                    _tracksWithNHitsContainer.getTracksWithNHitsVec( nHits ).push_back( trackAR );
                  }
                }
              }
            }
          }       
        }
      }
    }
  }
}


int SiliconTrackingAlg::BuildTrackFTD(TrackExtended * trackAR, int * nLR, int iS) {
  //  std::cout << "BuildTrackFTD: Layers = " << nLR[0] << " " << nLR[1] << " " << nLR[2] << endmsg;
  
  // initialise a helix from the track
  HelixClass helix;
  const float d0 = trackAR->getD0();
  const float z0 = trackAR->getZ0();
  const float phi0 = trackAR->getPhi();
  const float tanlambda = trackAR->getTanLambda();
  const float omega = trackAR->getOmega();
  helix.Initialize_Canonical(phi0,d0,z0,omega,tanlambda,_bField);
  float ref[3] = {helix.getReferencePoint()[0],
    helix.getReferencePoint()[1],
    helix.getReferencePoint()[2]};
  
  for (int iL=0; iL < static_cast<int>(_nlayersFTD); ++iL) {
    if (iL != nLR[0] && iL != nLR[1] && iL != nLR[2]) {
      float point[3];
      float ZL = _zLayerFTD[iL];
      if (iS == 0) 
        ZL = - ZL;
      helix.getPointInZ(ZL,ref,point);
      //      float Phi = atan2(point[1],point[0]);
      //      int iPhi = int(Phi/_dPhiFTD);
      float distMin = 1e+6;
      TrackerHitExtended * attachedHit = NULL;
      for (int ip=0;ip<=_nDivisionsInPhiFTD;++ip) {
        int iP = ip;
        if (iP < 0)
          iP = ip + _nDivisionsInPhiFTD;
        if (iP >= _nDivisionsInPhiFTD)
          iP = ip - _nDivisionsInPhiFTD;   
        int iCode = iS + 2*iL + 2*_nlayersFTD*iP;
        TrackerHitExtendedVec& hitVec = _sectorsFTD[iCode];
        int nH = int(hitVec.size());
        for (int iH=0; iH<nH; ++iH) {
          TrackerHitExtended * hit = hitVec[iH];
	  edm4hep::ConstTrackerHit trkHit = hit->getTrackerHit();
          float pos[3];
          for (int i=0;i<3;++i)
            pos[i] = float(trkHit.getPosition()[i]);
          float distance[3];
          float time = helix.getDistanceToPoint(pos,distance);
          if (time < 1.0e+10) {
            if (distance[2] < distMin) {
              distMin = distance[2];
              attachedHit = hit;
            }
          }
        }
      }
      //      std::cout << "Layer = " << iL << "  distMin = " << distMin << endmsg;
      if (distMin < _minDistCutAttach && attachedHit != NULL) {
        int iopt = 2;
        AttachHitToTrack( trackAR, attachedHit, iopt);
      }
    }
  }
  TrackerHitExtendedVec& hitVec = trackAR->getTrackerHitExtendedVec();
  int nH = int (hitVec.size());
  return nH;
}

int SiliconTrackingAlg::AttachHitToTrack(TrackExtended * trackAR, TrackerHitExtended * hit, int iopt) {
  
  int attached = 0;
  TrackerHitExtendedVec& hitVec = trackAR->getTrackerHitExtendedVec();
  int nHits = int(hitVec.size());
  
  double * xh = new double[nHits+1];
  double * yh = new double[nHits+1];
  float  * zh = new float[nHits+1];
  double * wrh = new double[nHits+1];
  float * wzh = new float[nHits+1];
  float * rh = new float[nHits+1];
  float * ph = new float[nHits+1];
  float par[5];
  float epar[15];
  
  for (int i=0; i<nHits; ++i) {
    edm4hep::ConstTrackerHit trkHit = hitVec[i]->getTrackerHit();
    xh[i] = double(trkHit.getPosition()[0]);
    yh[i] = double(trkHit.getPosition()[1]);
    zh[i] = float(trkHit.getPosition()[2]);
    ph[i] = float(atan2(yh[i],xh[i]));
    rh[i] = float(sqrt(xh[i]*xh[i]+yh[i]*yh[i]));
    float rR = hitVec[i]->getResolutionRPhi();
    float rZ = hitVec[i]->getResolutionZ();
    wrh[i] = double(1.0/(rR*rR));
    wzh[i] = 1.0/(rZ*rZ);
  }
  
  edm4hep::ConstTrackerHit trkHit = hit->getTrackerHit();
  xh[nHits] = double(trkHit.getPosition()[0]);
  yh[nHits] = double(trkHit.getPosition()[1]);
  zh[nHits] = float(trkHit.getPosition()[2]);
  ph[nHits] = float(atan2(yh[nHits],xh[nHits]));
  rh[nHits] = float(sqrt(xh[nHits]*xh[nHits]+yh[nHits]*yh[nHits]));
  
  float rR = hit->getResolutionRPhi();
  float rZ = hit->getResolutionZ();
  wrh[nHits] = double(1.0/(rR*rR));
  wzh[nHits] = 1.0/(rZ*rZ);
  
  
  int NPT = nHits + 1;
  
  // SJA:FIXME the newtonian part is giving crazy results for FTD so just use iopt 2 for simply attaching hits 
  // using SIT and VXD doesn't seem to give any problems, so make it a function parameter and let the caller decide
  //  int iopt = 3;
  
  float chi2RPhi = 0 ;
  float chi2Z = 0 ;
  
  
  int error = _fastfitter->fastHelixFit(NPT, xh, yh, rh, ph, wrh, zh, wzh,iopt, par, epar, chi2RPhi, chi2Z);
  par[3] = par[3]*par[0]/fabs(par[0]);
  
  
  float omega = par[0];
  float tanlambda = par[1];
  float phi0 = par[2];
  float d0 = par[3];
  float z0 = par[4];
  float chi2 = FLT_MAX;
  int ndf = INT_MAX;
  
  if (NPT == 3) {
    chi2 = chi2RPhi*_chi2WRPhiTriplet+chi2Z*_chi2WZTriplet;
  }
  if (NPT == 4) {
    chi2 = chi2RPhi*_chi2WRPhiQuartet+chi2Z*_chi2WZQuartet;
  }
  if (NPT > 4) {
    chi2 = chi2RPhi*_chi2WRPhiSeptet+chi2Z*_chi2WZSeptet;
  }
  ndf = 2*NPT-5;
  
  
  if ( error == 0 && chi2/float(ndf) < _chi2FitCut ) {
    trackAR->addTrackerHitExtended(hit);
    hit->addTrackExtended( trackAR );
    trackAR->setChi2( chi2 );
    trackAR->setOmega( omega );
    trackAR->setTanLambda( tanlambda );
    trackAR->setD0( d0 );
    trackAR->setZ0( z0 );
    trackAR->setPhi( phi0 );
    trackAR->setNDF( ndf );
    trackAR->setCovMatrix( epar );
    attached = 1;
    debug() << "Attachement succeeded chi2/float(ndf) = " << chi2/float(ndf) << "  cut = " <<  _chi2FitCut  << " chi2RPhi = " << chi2RPhi << " chi2Z = " << chi2Z << " error = " << error << endmsg;
  } else {
    debug() << "Attachement failed chi2/float(ndf) = " << chi2/float(ndf) << "  cut = " <<  _chi2FitCut  << " chi2RPhi = " << chi2RPhi << " chi2Z = " << chi2Z << " error = " << error << endmsg;
  }
  
  delete[] xh;
  delete[] yh;
  delete[] zh;
  delete[] wrh;
  delete[] wzh;
  delete[] rh;
  delete[] ph;
  
  return attached;
  
  
}

void SiliconTrackingAlg::FinalRefit(edm4hep::TrackCollection* trk_col) {
  
  int nTracks = int(_trackImplVec.size());
  
  int nSiSegments = 0;        
  float eTot = 0.;
  float pxTot = 0.;
  float pyTot = 0.;
  float pzTot = 0.;
  debug() << "Total " << nTracks << " candidate tracks will be dealed" << endmsg;
  for (int iTrk=0;iTrk<nTracks;++iTrk) {
    
    TrackExtended * trackAR = _trackImplVec[iTrk];    
    TrackerHitExtendedVec& hitVec = trackAR->getTrackerHitExtendedVec();
    
    int nHits = int(hitVec.size());
    //debug() << "  Track " << iTrk << " has hits: " << nHits << endmsg;
    if( nHits >= _minimalHits) {
      //    int * lh = new int[nHits];
      std::vector<int> lh;
      lh.resize(nHits);
      
      for (int i=0; i<nHits; ++i) {
        lh[i]=0;
      }
      
      float d0 = trackAR->getD0();
      float z0 = trackAR->getZ0();
      float omega = trackAR->getOmega();
      float tanlambda = trackAR->getTanLambda();
      float phi0 = trackAR->getPhi();
      
      HelixClass * helix = new HelixClass();
      helix->Initialize_Canonical(phi0, d0, z0, omega, 
                                  tanlambda, _bField);
      
      // get the point of closest approach to the reference point
      // here it is implicitly assumed that the reference point is the origin 
      float Pos[3];
      Pos[0] = -d0*sin(phi0);
      Pos[1] = d0*cos(phi0);
      Pos[2] = z0;
      
      // at this point is is possible to have hits from the same layer ...
      // so a check is made to ensure that the hit with the smallest distance to the 
      // current helix hypothosis is used, the other hit has lh set to 0 
      
      // start loop over the hits to
      for (int ihit=0;ihit<nHits;++ihit) {
        
        lh[ihit] = 1; // only hits which have lh=1 will be used for the fit
        
        // get the pointer to the lcio trackerhit for this hit
	edm4hep::ConstTrackerHit trkHit = hitVec[ihit]->getTrackerHit();
        
        int det = getDetectorID(trkHit);
        
        if (det == lcio::ILDDetID::VXD || det == lcio::ILDDetID::FTD || det == lcio::ILDDetID::SIT) { // only accept VXD, FTD or SIT
          
          
          //        int layer = getLayerID(trkHit);
          //        int moduleIndex = getModuleID(trkHit);
          
          // start a double loop over the hits which have already been checked 
          for (int lhit=0;lhit<ihit;++lhit) {
            
            // get the pointer to the lcio trackerhit for the previously checked hit
	    edm4hep::ConstTrackerHit trkHitS = hitVec[lhit]->getTrackerHit();
            
            
            //          int layerS = getLayerID(trkHitS);
            //          int moduleIndexS = getModuleID(trkHitS);
            
            // SJA:FIXME: check to see if allowing no hits in the same sensor vs no hits in the same layer works 
            // if they are on the same layer and the previously checked hits has been declared good for fitting
            //          if ((trkHitS->getType() == trkHit->getType()) && (lh[lhit] == 1)) {
            // check if the hits have the same layer and petal number
            //          hitVec[ihit]->
            //          if ((layer == layerS) && (moduleIndex==moduleIndexS) && (lh[lhit] == 1)) {
            if ( (trkHit.getCellID() == trkHitS.getCellID()) && (lh[lhit] == 1)) {
              
              // get the position of the hits 
              float xP[3];
              float xPS[3];
              for (int iC=0;iC<3;++iC) {
                xP[iC] = float(trkHit.getPosition()[iC]);
                xPS[iC] = float(trkHitS.getPosition()[iC]);
              }
              
              // get the intersection of the helix with the either the cylinder or plane containing the hit
              float Point[6];
              float PointS[6];
              
              if (det == lcio::ILDDetID::FTD) {

                float time = helix->getPointInZ(xP[2],Pos,Point);
                time = helix->getPointInZ(xPS[2],Pos,PointS);

              } else {

                float RAD = sqrt(xP[0]*xP[0]+xP[1]*xP[1]);
                float RADS = sqrt(xPS[0]*xPS[0]+xPS[1]*xPS[1]);
                float time = helix->getPointOnCircle(RAD,Pos,Point);
                time = helix->getPointOnCircle(RADS,Pos,PointS);

              }
              
              float DIST = 0;
              float DISTS = 0;
              
              // get the euclidean distance between the hit and the point of intersection
              for (int iC=0;iC<3;++iC) {
                DIST += (Point[iC]-xP[iC])*(Point[iC]-xP[iC]);
                DISTS += (PointS[iC]-xPS[iC])*(PointS[iC]-xPS[iC]);
              }
              if (DIST < DISTS) {
                lh[lhit] = 0;
              }
              else {
                lh[ihit] = 0;
              }
              break;
            }
          }
        }
      }
      
      delete helix;
      
      std::vector<ConstTrackerHit> trkHits;
      std::vector<ConstTrackerHit> trkHits_used_inFit;
      
      int nFit = 0;
      for (int i=0; i<nHits; ++i) {
        // check if the hit has been rejected as being on the same layer and further from the helix lh==0
        if (lh[i] == 1) {
	  edm4hep::ConstTrackerHit trkHit = hitVec[i]->getTrackerHit();
	  debug() << "TrackerHit " << i << " id = " << trkHit.id() << endmsg;
          nFit++;
          if(trkHit.isAvailable()) { 
            trkHits.push_back(trkHit);   
          }
          else{
            throw EVENT::Exception( std::string("SiliconTrackingAlg::FinalRefit: TrackerHit pointer == NULL ")  ) ;
          }
        }
        else { // reject hit 
               // SJA:FIXME missuse of type find a better way to signal rejected hits
          hitVec[i]->setType(int(0));
        }
      }
      
      if( trkHits.size() < 3 ) {
        debug() << "SiliconTrackingAlg::FinalRefit: Cannot fit less than 3 hits. Number of hits =  " << trkHits.size() << endmsg;
        continue ; 
      }
      //TrackImpl* Track = new TrackImpl ;
      //auto track = trk_col->create();
      //fucd
      edm4hep::Track track;// = new edm4hep::Track;
      // setup initial dummy covariance matrix
      //std::vector<float> covMatrix;
      //covMatrix.resize(15);
      std::array<float,15> covMatrix;

      for (unsigned icov = 0; icov<covMatrix.size(); ++icov) {
        covMatrix[icov] = 0;
      }
      
      covMatrix[0]  = ( _initialTrackError_d0    ); //sigma_d0^2
      covMatrix[2]  = ( _initialTrackError_phi0  ); //sigma_phi0^2
      covMatrix[5]  = ( _initialTrackError_omega ); //sigma_omega^2
      covMatrix[9]  = ( _initialTrackError_z0    ); //sigma_z0^2
      covMatrix[14] = ( _initialTrackError_tanL  ); //sigma_tanl^2
      
      
      std::vector< std::pair<float, edm4hep::ConstTrackerHit> > r2_values;
      r2_values.reserve(trkHits.size());
      
      for (std::vector<edm4hep::ConstTrackerHit>::iterator it=trkHits.begin(); it!=trkHits.end(); ++it) {
        edm4hep::ConstTrackerHit h = *it;
        float r2 = h.getPosition()[0]*h.getPosition()[0]+h.getPosition()[1]*h.getPosition()[1];
        r2_values.push_back(std::make_pair(r2, *it));
      }
      
      sort(r2_values.begin(),r2_values.end());
      
      trkHits.clear();
      trkHits.reserve(r2_values.size());

      for (std::vector< std::pair<float, edm4hep::ConstTrackerHit> >::iterator it=r2_values.begin(); it!=r2_values.end(); ++it) {
        trkHits.push_back(it->second);
      }
      //std::cout << "fucd------------------3 " << _trksystem << std::endl;
      //for (unsigned ihit_indx=0 ; ihit_indx < trkHits.size(); ++ihit_indx) {
      //  std::cout << "fucd trk hit " << *trkHits[ihit_indx] << " " << trkHits[ihit_indx]->getCovMatrix()[0]
      //		  << " " << BitSet32(trkHits[ihit_indx]->getType())[ ILDTrkHitTypeBit::COMPOSITE_SPACEPOINT ] << endmsg;
      //}
      /*
      auto _trackSystemSvc = service<ITrackSystemSvc>("TrackSystemSvc");
      if ( !_trackSystemSvc ) {
	error() << "Failed to find TrackSystemSvc ..." << endmsg;
	return;
      }
      _trksystem =  _trackSystemSvc->getTrackSystem();

      if( _trksystem == 0 ){
	error() << "Cannot initialize MarlinTrkSystem of Type: KalTest" <<endmsg;
	return;
      }
      debug() << "_trksystem pointer " << _trksystem << endmsg;
      
      _trksystem->setOption( IMarlinTrkSystem::CFG::useQMS,        _MSOn ) ;
      _trksystem->setOption( IMarlinTrkSystem::CFG::usedEdx,       _ElossOn) ;
      _trksystem->setOption( IMarlinTrkSystem::CFG::useSmoothing,  _SmoothOn) ;
      _trksystem->init() ;
      */
      bool fit_backwards = IMarlinTrack::backward;
      
      MarlinTrk::IMarlinTrack* marlinTrk = nullptr;
      try{
	marlinTrk = _trksystem->createTrack();
      }
      catch(...){
	error() << "Cannot create MarlinTrack ! " << endmsg;
	return;
      }
      
      int status = 0;
      debug() << "call createFinalisedLCIOTrack now" << endmsg;
      try {
        status = MarlinTrk::createFinalisedLCIOTrack(marlinTrk, trkHits, &track, fit_backwards, covMatrix, _bField, _maxChi2PerHit);
      } catch (...) {
	//      delete Track;
        //      delete marlinTrk;
        error() << "MarlinTrk::createFinalisedLCIOTrack fail" << endmsg;
        //throw ;
      }
      debug() << "createFinalisedLCIOTrack finish" << endmsg;
      /*
#ifdef MARLINTRK_DIAGNOSTICS_ON
      if ( status != IMarlinTrack::success && _runMarlinTrkDiagnostics ) {        
        void * dcv = _trksystem->getDiagnositicsPointer();
        DiagnosticsController* dc = static_cast<DiagnosticsController*>(dcv);
        dc->skip_current_track();
      }        
#endif
      */
      
      std::vector<std::pair<edm4hep::ConstTrackerHit , double> > hits_in_fit ;  
      std::vector<std::pair<edm4hep::ConstTrackerHit , double> > outliers ;
      std::vector<edm4hep::ConstTrackerHit> all_hits;    
      all_hits.reserve(300);
      
      marlinTrk->getHitsInFit(hits_in_fit);
      
      for ( unsigned ihit = 0; ihit < hits_in_fit.size(); ++ihit) {
	debug() << "Hit id =" << hits_in_fit[ihit].first.id() << endmsg;
	edm4hep::ConstTrackerHit trk = hits_in_fit[ihit].first;
        all_hits.push_back(trk);//hits_in_fit[ihit].first);
      }
      
      UTIL::BitField64 cellID_encoder( UTIL::ILDCellID0::encoder_string ) ; 
      
      MarlinTrk::addHitNumbersToTrack(&track, all_hits, true, cellID_encoder);
      
      marlinTrk->getOutliers(outliers);
      
      for ( unsigned ihit = 0; ihit < outliers.size(); ++ihit) {
        all_hits.push_back(outliers[ihit].first);
      }
      
      MarlinTrk::addHitNumbersToTrack(&track, all_hits, false, cellID_encoder);
      
      delete marlinTrk;

      int nhits_in_vxd = track.getSubDetectorHitNumbers(0);
      int nhits_in_ftd = track.getSubDetectorHitNumbers(1);
      int nhits_in_sit = track.getSubDetectorHitNumbers(2);
      
      //debug() << " Hit numbers for Track "<< track.id() << ": "
      debug() << " Hit numbers for Track "<< iTrk <<": "
	      << " vxd hits = " << nhits_in_vxd
	      << " ftd hits = " << nhits_in_ftd
	      << " sit hits = " << nhits_in_sit
	      << endmsg;
      
      //if (nhits_in_vxd > 0) Track->setTypeBit( lcio::ILDDetID::VXD ) ;
      //if (nhits_in_ftd > 0) Track->setTypeBit( lcio::ILDDetID::FTD ) ;
      //if (nhits_in_sit > 0) Track->setTypeBit( lcio::ILDDetID::SIT ) ;

      if( status != IMarlinTrack::success ) {       
        //delete track;
        debug() << "FinalRefit: Track fit failed with error code " << status << " track dropped. Number of hits = "<< trkHits.size() << endmsg;       
        continue ;
      }
      
      if( track.getNdf() < 0) {       
        //delete track;
        debug() << "FinalRefit: Track fit returns " << track.getNdf() << " degress of freedom track dropped. Number of hits = "<< trkHits.size() << endmsg;       
	//delete track;
        continue ;
      }
      
      //trk_col->addElement(Track);     
      //fucd
      trk_col->push_back(track);
      for(int i=0;i<track.trackStates_size();i++){
	// 1 = lcio::EVENT::TrackState::AtIP
	edm4hep::TrackState trkStateIP = track.getTrackStates(i);
	if(trkStateIP.location !=1) continue;
      /*
      if (trkStateIP == 0) {
        debug() << "SiliconTrackingAlg::FinalRefit: Track fit returns " << track->getNdf() << " degress of freedom track dropped. Number of hits = "<< trkHits.size() << endmsg;       
        throw EVENT::Exception( std::string("SiliconTrackingAlg::FinalRefit: trkStateIP pointer == NULL ")  ) ;
      }
      */
      // note trackAR which is of type TrackExtended, only takes fits set for ref point = 0,0,0
	trackAR->setOmega(trkStateIP.omega);
	trackAR->setTanLambda(trkStateIP.tanLambda);
	trackAR->setPhi(trkStateIP.phi);
	trackAR->setD0(trkStateIP.D0);
	trackAR->setZ0(trkStateIP.Z0);
	
	float cov[15];
	
	for (int i = 0 ; i<15 ; ++i) {
	  cov[i] = trkStateIP.covMatrix.operator[](i);
	}
      
	trackAR->setCovMatrix(cov);
	trackAR->setChi2(track.getChi2());
	trackAR->setNDF(track.getNdf());
      
	nSiSegments++;
      
	HelixClass helix_final;
      
	helix_final.Initialize_Canonical(trkStateIP.phi,trkStateIP.D0,trkStateIP.Z0,trkStateIP.omega,trkStateIP.tanLambda,_bField);
      
	float trkPx = helix_final.getMomentum()[0];
	float trkPy = helix_final.getMomentum()[1];
	float trkPz = helix_final.getMomentum()[2];
	float trkP = sqrt(trkPx*trkPx+trkPy*trkPy+trkPz*trkPz);
	eTot += trkP;
	pxTot += trkPx;
	pyTot += trkPy;
	pzTot += trkPz;
      }
    }
  }
  
  debug() << " -> run " << _nRun << " event " << _nEvt << endmsg;
  debug() << "Number of Si tracks = " << nSiSegments << endmsg;
  debug() << "Total 4-momentum of Si Tracks : E = " << std::setprecision(7) << eTot
	  << " Px = " << pxTot
	  << " Py = " << pyTot
	  << " Pz = " << pzTot << endmsg;
}

StatusCode SiliconTrackingAlg::setupGearGeom(){
  auto _gear = service<IGearSvc>("GearSvc");
  if ( !_gear ) {
    error() << "Failed to find GearSvc ..." << endmsg;
    return StatusCode::FAILURE;
  }
  gear::GearMgr* gearMgr = _gear->getGearMgr();
  _bField = gearMgr->getBField().at( gear::Vector3D( 0.,0.,0.)  ).z() ;
  debug() << "Field " << _bField << endmsg;
  //-- VXD Parameters--
  _nLayersVTX = 0 ;
  const gear::VXDParameters* pVXDDetMain = 0;
  const gear::VXDLayerLayout* pVXDLayerLayout = 0;
  
  try{
    
    debug() << " filling VXD parameters from gear::SITParameters " << endmsg ;
    
    pVXDDetMain = &gearMgr->getVXDParameters();
    pVXDLayerLayout = &(pVXDDetMain->getVXDLayerLayout());
    _nLayersVTX = pVXDLayerLayout->getNLayers();
  }
  catch( gear::UnknownParameterException& e){
    
    debug() << " ### gear::VXDParameters Not Present in GEAR FILE" << endmsg ;
    
  }
  
  
  
  //-- SIT Parameters--
  _nLayersSIT = 0 ;
  const gear::ZPlanarParameters* pSITDetMain = 0;
  const gear::ZPlanarLayerLayout* pSITLayerLayout = 0;
  
  try{
    
    debug() << " filling SIT parameters from gear::SITParameters " << endmsg ;
    
    pSITDetMain = &gearMgr->getSITParameters();
    pSITLayerLayout = &(pSITDetMain->getZPlanarLayerLayout());
    _nLayersSIT = pSITLayerLayout->getNLayers();
    
  }
  catch( gear::UnknownParameterException& e){
    
    debug() << " ### gear::SITParameters Not Present in GEAR FILE" << endmsg ;
    
  }
  
  if( _nLayersSIT == 0 ){
    // try the old LOI style key value pairs as defined in the SSit03 Mokka drive
    try{
      
      info() << "  SiliconTrackingAlg - Simple Cylinder Based SIT using parameters defined by SSit03 Mokka driver " << endmsg ;
      
      // SIT
      
      const gear::GearParameters& pSIT = gearMgr->getGearParameters("SIT");
      
      const std::vector<double>& SIT_r   =  pSIT.getDoubleVals("SITLayerRadius" )  ;
      const std::vector<double>& SIT_hl  =  pSIT.getDoubleVals("SITSupportLayerHalfLength" )  ;
      
      _nLayersSIT = SIT_r.size() ; 
      
      if (_nLayersSIT != SIT_r.size() || _nLayersSIT != SIT_hl.size()) {
        
        error() << "ILDSITCylinderKalDetector miss-match between DoubleVec and nlayers exit(1) called from file " << __FILE__ << " line " << __LINE__  << endmsg ;
        exit(1);
        
      }
    }
    catch( gear::UnknownParameterException& e){
      
      debug() << " ### gear::SIT Parameters from as defined in SSit03 Not Present in GEAR FILE" << endmsg ;
      
    } 
    
  }
  
  
  
  //-- FTD Parameters--
  _petalBasedFTDWithOverlaps = false;  
  _nlayersFTD = 0;
  
  try{
    
    debug() << " filling FTD parameters from gear::FTDParameters " << endmsg ;
    
    const gear::FTDParameters&   pFTD      = gearMgr->getFTDParameters();
    const gear::FTDLayerLayout&  ftdlayers = pFTD.getFTDLayerLayout() ;
    
    _nlayersFTD = ftdlayers.getNLayers() ;
    
    for (unsigned int disk=0; disk < _nlayersFTD; ++disk) {
      
      _zLayerFTD.push_back( ftdlayers.getSensitiveZposition(disk, 0, 1) ); // front petal even numbered
      
      if ( ftdlayers.getNPetals(disk) > 0) {
        _zLayerFTD.push_back( ftdlayers.getSensitiveZposition(disk, 1, 1) );  // front petal odd numbered
        _petalBasedFTDWithOverlaps = true;
      }
      
    }
    
    // SJA: Here we increase the size of _nlayersFTD as we are treating the 
    _nlayersFTD =_zLayerFTD.size() ;     
    
  }
  catch( gear::UnknownParameterException& e){
    
    debug() << " ### gear::FTDParameters Not Present in GEAR FILE" << endmsg ;
    
  } 
  
  if( _nlayersFTD == 0 ){
    
    // FTD
    try{
      
      info() << "  SiliconTrackingAlg - Simple Disc Based FTD using parameters defined by SFtd05 Mokka driver " << endmsg ;
      
      const gear::GearParameters& pFTD = gearMgr->getGearParameters("FTD");
      
      const std::vector<double>* pFTD_z   = NULL;
      
      info() << " For FTD using parameters defined by SFtd05 Mokka driver " << endmsg ;
      
      pFTD_z = &pFTD.getDoubleVals("FTDZCoordinate" )  ;
      
      _nlayersFTD = pFTD_z->size();
      
      for (unsigned int i = 0; i<_nlayersFTD; ++i) {
        _zLayerFTD.push_back((*pFTD_z)[i]);
      }
    }
    catch( gear::UnknownParameterException& e){
      
      debug() << " ### gear::FTD Parameters as defined in SFtd05 Not Present in GEAR FILE" << endmsg ;
      
    } 
  }
  return StatusCode::SUCCESS;
}

void SiliconTrackingAlg::TracksWithNHitsContainer::clear()
{
  for (std::vector< TrackExtendedVec >::iterator trackVecIter = _tracksNHits.begin();
       trackVecIter < _tracksNHits.end(); trackVecIter++)
  {
    for (TrackExtendedVec::iterator trackIter = trackVecIter->begin();
         trackIter < trackVecIter->end(); trackIter++)
    {
      delete *trackIter;
    }
    
    trackVecIter->clear();
  }
}

