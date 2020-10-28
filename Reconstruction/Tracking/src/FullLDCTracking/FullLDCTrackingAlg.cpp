#include "FullLDCTrackingAlg.h"

#include "DataHelper/Navigation.h"
#include "Tracking/TrackingHelper.h"

#include <GearSvc/IGearSvc.h>

#include <edm4hep/TrackerHitConst.h>
#include <edm4hep/TrackerHit.h>
#include <edm4hep/TrackConst.h>
#include <edm4hep/Track.h>

#include <iostream>
#include <algorithm>
#include <memory>

#include <math.h>
#include <map>

//#include "DataHelper/ClusterShapes.h"
#include <gear/GEAR.h>
#include <gear/GearParameters.h>
#include <gear/BField.h>
#include <gear/VXDLayerLayout.h>
#include <gear/VXDParameters.h>
#include "gear/FTDLayerLayout.h"
#include "gear/FTDParameters.h"
#include <gear/TPCParameters.h>
#include <gear/PadRowLayout2D.h>

#include "TrackSystemSvc/IMarlinTrack.h"
#include "TrackSystemSvc/IMarlinTrkSystem.h"
#include "TrackSystemSvc/MarlinTrkUtils.h"
#include "TrackSystemSvc/LCIOTrackPropagators.h"

#include <UTIL/LCTOOLS.h>
//#include <UTIL/LCRelationNavigator.h>

// #include "MarlinTrk/MarlinTrkDiagnostics.h"
#ifdef MARLINTRK_DIAGNOSTICS_ON
#include "MarlinTrk/DiagnosticsController.h"
#endif

#include <UTIL/BitField64.h>
//#include <UTIL/BitSet32.h>
#include <UTIL/ILDConf.h>

#include <climits>
#include <cmath>

#include "gsl/gsl_randist.h"
#include "gsl/gsl_cdf.h"

#include <vector>
#include <bitset>

typedef std::vector<edm4hep::ConstTrackerHit> ConstTrackerHitVec;

using namespace edm4hep ;
using namespace MarlinTrk ;

/** debug printout helper method */
std::string toString( int iTrk, edm4hep::ConstTrack tpcTrack, float bField=3.5 ) {
  
  int   nHits    = int( tpcTrack.trackerHits_size() );
  float d0TPC    = getD0(tpcTrack);
  float z0TPC    = getZ0(tpcTrack);
  float omegaTPC = getOmega(tpcTrack);
  float phi0TPC  = getPhi(tpcTrack);
  float tanLTPC  = getTanLambda(tpcTrack);
  float Chi2TPC  = tpcTrack.getChi2()/float(tpcTrack.getNdf());
  int   ndfTPC   = tpcTrack.getNdf();

  int nlinkedTracks = tpcTrack.tracks_size();
  
  
  HelixClass helixTPC;
  
  helixTPC.Initialize_Canonical(phi0TPC,d0TPC,z0TPC,omegaTPC,tanLTPC, bField);
  
  char strg[200];
  
  float pxTPC = helixTPC.getMomentum()[0];
  float pyTPC = helixTPC.getMomentum()[1];
  float pzTPC = helixTPC.getMomentum()[2];
  const float ptot = sqrt(pxTPC*pxTPC+pyTPC*pyTPC+pzTPC*pzTPC);

  sprintf(strg,"%3i   %5i %9.3f  %9.3f  %9.3f  %7.2f  %7.2f  %7.2f %4i %4i %8.3f %8i",iTrk,tpcTrack.id(),
	  ptot, d0TPC,z0TPC,pxTPC,pyTPC,pzTPC,nHits,ndfTPC,Chi2TPC,nlinkedTracks);

  return std::string( strg ) ;
}

DECLARE_COMPONENT( FullLDCTrackingAlg )
FullLDCTrackingAlg::FullLDCTrackingAlg(const std::string& name, ISvcLocator* svcLoc) : GaudiAlgorithm(name, svcLoc) {  
  // _description = "Performs full tracking in ILD detector" ;  
  
  _encoder = new UTIL::BitField64(lcio::ILDCellID0::encoder_string);
  
  // Input tracker hit collections
  declareProperty("FTDPixelTrackerHits", _FTDPixelTrackerHitColHdl, "handler of FTD Pixel Hit Collection");
  declareProperty("FTDSpacePoints", _FTDSpacePointColHdl, "FTD FTDSpacePoint Collection");
  declareProperty("VTXTrackerHits", _VTXTrackerHitColHdl, "VTX Hit Collection");
  declareProperty("SITTrackerHits", _SITTrackerHitColHdl, "SIT Hit Collection");
  declareProperty("SETTrackerHits", _SETTrackerHitColHdl, "SET Hit Collection");
  //declareProperty("ETDTrackerHits", _ETDTrackerHitColHdl, "ETD Hit Collection");
  declareProperty("TPCTrackerHits", _TPCTrackerHitColHdl, "TPC Hit Collection");
  declareProperty("SITRawHits", _inSITRawColHdl, "SIT Raw Hit Collection of SpacePoints");
  declareProperty("SETRawHits", _inSETRawColHdl, "SET Raw Hit Collection of SpacePoints");
  declareProperty("FTDRawHits", _inFTDRawColHdl, "FTD Raw Hit Collection of SpacePoints");
  //declareProperty("VTXRawHits", _inVXDRawColHdl, "VXD SimTrackerHit collection");
  
  // Input track collections
  declareProperty("TPCTracks", _TPCTrackColHdl, "TPC Track Collection"); 
  declareProperty("SiTracks", _SiTrackColHdl, "Si Track Collection");
  
  // Input relation collections
  
  /*
  registerInputCollection(LCIO::LCRELATION,
                          "TPCTracksMCPRelColl",
                          "TPC Track to MCP Relation Collection Name",
                          _TPCTrackMCPCollName,
                          std::string("TPCTracksMCP"));
  
  registerInputCollection(LCIO::LCRELATION,
                          "SiTracksMCPRelColl",
                          "Si Track to Collection",
                          _SiTrackMCPCollName,
                          std::string("SiTracksMCP"));
  */ 
  // Output track collection
  declareProperty("OutputTracks", _OutputTrackColHdl, "Full LDC track collection name");
}



StatusCode FullLDCTrackingAlg::initialize() { 
  
  // usually a good idea to
  // printParameters();  
  _nRun = -1 ;
  _nEvt = 0 ;
  PI = acos(-1.);
  PIOVER2 = 0.5*PI;
  TWOPI = 2*PI;
  
  // set up the geometery needed by KalTest
  //FIXME: for now do KalTest only - make this a steering parameter to use other fitters
  auto _trackSystemSvc = service<ITrackSystemSvc>("TrackSystemSvc");
  if ( !_trackSystemSvc ) {
    error() << "Fail to find TrackSystemSvc ..." << endmsg;
  } 
  
  _trksystem = _trackSystemSvc->getTrackSystem(this);
  if( _trksystem == 0 ){
    error() << "Cannot initialize MarlinTrkSystem of Type: KalTest" <<endmsg;
    return StatusCode::FAILURE;
  }

  _trksystem->setOption( MarlinTrk::IMarlinTrkSystem::CFG::useQMS,        _MSOn ) ;
  _trksystem->setOption( MarlinTrk::IMarlinTrkSystem::CFG::usedEdx,       _ElossOn) ;
  _trksystem->setOption( MarlinTrk::IMarlinTrkSystem::CFG::useSmoothing,  _SmoothOn) ;
  _trksystem->init() ;
    
  this->setupGearGeom();
  return GaudiAlgorithm::initialize();
}

/*
void FullLDCTrackingAlg::processRunHeader( LCRunHeader* run) { 
  
  _nRun++ ;
  _nEvt = 0;
  streamlog_out(DEBUG5) << endmsg;
  streamlog_out(DEBUG5) << "FullLDCTrackingAlg ---> new run : run number = " << run->getRunNumber() << endmsg;
  
} 
*/

StatusCode FullLDCTrackingAlg::execute() { 
  
  // debug() << endmsg;
  debug() << "FullLDCTrackingAlg -> run = " << 0/*evt->getRunNumber()*/ << "  event = " << _nEvt << endmsg;
  // debug() << endmsg;
  auto outCol = _OutputTrackColHdl.createAndPut();
  
  prepareVectors();
  debug() << "************************************PrepareVectors done..." << endmsg;

  debug() << "************************************Merge TPC/Si ..." << endmsg;
  
  MergeTPCandSiTracks();
  debug() << "************************************Merging done ..." << endmsg;

  MergeTPCandSiTracksII();
  debug() << "************************************Merging II done ..." << endmsg;

  Sorting(_allCombinedTracks);
  debug() << "************************************Sorting by Chi2/NDF done ..." << endmsg;

  debug() << "************************************Selection of all 2 track combininations ..." << endmsg;
  SelectCombinedTracks();
  debug() << "************************************Selection of all 2 track combininations done ..." << endmsg;

  debug() << "************************************Trying non combined tracks ..." << endmsg;
  AddNotCombinedTracks( );
  debug() << "************************************Non combined tracks added ..." << endmsg;
  CheckTracks( );

  debug() << "************************************Add Non assigned hits ..." << endmsg;
  AddNotAssignedHits();
  debug() << "************************************Non assigned hits added ..." << endmsg;

  AddTrackColToEvt(_trkImplVec, outCol);
  debug() << "Collections added to event ..." << endmsg;
  CleanUp();
  debug() << "Cleanup is done." << endmsg;
  _nEvt++;
  //  getchar();
  // streamlog_out(DEBUG5) << endmsg;
  // streamlog_out(DEBUG5) << endmsg;
  return StatusCode::SUCCESS;
  
}

void FullLDCTrackingAlg::AddTrackColToEvt(TrackExtendedVec & trkVec, edm4hep::TrackCollection* colTRK) {
  
  //LCCollectionVec * colTRK = new LCCollectionVec(LCIO::TRACK);
  // if we want to point back to the hits we need to set the flag
  //LCFlagImpl trkFlag(0) ;
  //trkFlag.setBit( LCIO::TRBIT_HITS ) ;
  //colTRK->setFlag( trkFlag.getFlag()  ) ;  
    
  //  LCCollectionVec * colRel = NULL;
    
  int nTrkCand = int(trkVec.size());
  
  int nTotTracks = 0;
  float eTot  = 0.0;
  float pxTot = 0.0;
  float pyTot = 0.0;
  float pzTot = 0.0;
  
  //SJA:FIXME: So here we are going to do one final refit. This can certainly be optimised, but rather than worry about the mememory management right now lets make it work, and optimise it later ...
  
  
  for (int iTRK=0;iTRK<nTrkCand;++iTRK) {
    
    TrackExtended * trkCand = trkVec[iTRK];
    TrackerHitExtendedVec& hitVec = trkCand->getTrackerHitExtendedVec();
    
    std::vector<edm4hep::ConstTrackerHit> trkHits;
    
    debug() << " Trying to add track " << trkCand << " to final lcio collection " << endmsg;
        
    int nHits = int(hitVec.size());
    
    debug() << " Trying to add track " << trkCand << " to final lcio collection nHits = " << nHits << endmsg;
    
    for (int ihit=0;ihit<nHits;++ihit) {
    
      if( hitVec[ihit]->getUsedInFit() == false )  {
        debug() << "rejecting hit for track " << trkCand << " at zhit  " <<  hitVec[ihit]->getTrackerHit().getPosition()[2] << endmsg;
        continue;
      }

      edm4hep::ConstTrackerHit trkHit = hitVec[ihit]->getTrackerHit();
      
      if(trkHit.isAvailable()) {
        trkHits.push_back(trkHit);   
      }
      else{
        throw EVENT::Exception( std::string("FullLDCTrackingAlg::AddTrackColToEvt: TrackerHit pointer == NULL ")  ) ;
      }
      
    }
    
    
    if( trkHits.size() < 3 ) {
      debug() << "FullLDCTrackingAlg::AddTrackColToEvt: Cannot fit less than 3 hits. Number of hits =  " << trkHits.size() << endmsg;
      continue ; 
    }
    
    edm4hep::Track track;// = new edm4hep::Track;
    
    // setup initial dummy covariance matrix
    std::array<float,15> covMatrix;
    
    for (unsigned icov = 0; icov<covMatrix.size(); ++icov) {
      covMatrix[icov] = 0;
    }
    
    covMatrix[0]  = ( _initialTrackError_d0    ); //sigma_d0^2
    covMatrix[2]  = ( _initialTrackError_phi0  ); //sigma_phi0^2
    covMatrix[5]  = ( _initialTrackError_omega ); //sigma_omega^2
    covMatrix[9]  = ( _initialTrackError_z0    ); //sigma_z0^2
    covMatrix[14] = ( _initialTrackError_tanL  ); //sigma_tanl^2
    
    // get the track state at the last hit at the outer most hit
    
    GroupTracks * group = trkCand->getGroupTracks();

    edm4hep::TrackState ts_initial;
    
    bool prefit_set = false;
    
    debug() << "Track Group = " << group << endmsg;
    
    if( group ) debug() << "Track Group size = " << group->getTrackExtendedVec().size() << endmsg;

    if (group != NULL && group->getTrackExtendedVec().size() > 0) {
    
      // get the second track as this must be the one furthest from the IP

      TrackExtended* te = 0;
      
      if(group->getTrackExtendedVec().size()==1) {
        te = group->getTrackExtendedVec()[0];
      } else {
        te = group->getTrackExtendedVec()[1];
      }
      
      if(hasTrackStateAt(te->getTrack(), 3/*lcio::TrackState::AtLastHit*/)){

        debug() << "Initialise Fit with trackstate from last hit" << group << endmsg;
	
        ts_initial = getTrackStateAt(te->getTrack(), 3/*lcio::TrackState::AtLastHit*/);
                        
        prefit_set = true;
      }
    }
    
    if( !prefit_set ) { // use parameters at IP
      debug() << "Initialise Fit with trackstate from IP " << group << endmsg;
      
      ts_initial.D0        = trkCand->getD0();
      ts_initial.phi       = trkCand->getPhi();
      ts_initial.Z0        = trkCand->getZ0();
      ts_initial.omega     = trkCand->getOmega();
      ts_initial.tanLambda = trkCand->getTanLambda();
      
      edm4hep::Vector3f ref(0,0,0);
      
      ts_initial.referencePoint = ref;
      
      ts_initial.location = 1/*lcio::TrackState::AtIP*/;
    }
    
    ts_initial.covMatrix = covMatrix;
        
    // sort hits in R
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
    
    bool fit_backwards = IMarlinTrack::backward;
    
    MarlinTrk::IMarlinTrack* marlinTrk = _trksystem->createTrack();
    
    
    int error = 0;
    
    try {
      
      error = MarlinTrk::createFinalisedLCIOTrack(marlinTrk, trkHits, &track, fit_backwards, &ts_initial, _bField, _maxChi2PerHit);                              
      
    } catch (...) {
      
      //      delete track;
      //      delete marlinTrk;
      
      throw ;
      
    }
    
    
#ifdef MARLINTRK_DIAGNOSTICS_ON
    if ( error != IMarlinTrack::success && _runMarlinTrkDiagnostics ) {        
      void * dcv = _trksystem->getDiagnositicsPointer();
      DiagnosticsController* dc = static_cast<DiagnosticsController*>(dcv);
      dc->skip_current_track();
    }        
#endif
    
    
    std::vector<std::pair<edm4hep::ConstTrackerHit , double> > hits_in_fit ;  
    std::vector<std::pair<edm4hep::ConstTrackerHit , double> > outliers ;
    std::vector<edm4hep::ConstTrackerHit> all_hits;    
    all_hits.reserve(300);
    
    marlinTrk->getHitsInFit(hits_in_fit);
    
    for ( unsigned ihit = 0; ihit < hits_in_fit.size(); ++ihit) {
      all_hits.push_back(hits_in_fit[ihit].first);
    }
    
    UTIL::BitField64 cellID_encoder( lcio::ILDCellID0::encoder_string ) ; 
    
    MarlinTrk::addHitNumbersToTrack(&track, all_hits, true, cellID_encoder);
    
    marlinTrk->getOutliers(outliers);
    
    for ( unsigned ihit = 0; ihit < outliers.size(); ++ihit) {
      all_hits.push_back(outliers[ihit].first);
    }
    
    MarlinTrk::addHitNumbersToTrack(&track, all_hits, false, cellID_encoder);
    
    
    delete marlinTrk;
    
    if( error != IMarlinTrack::success ) {       
      debug() << "FullLDCTrackingAlg::AddTrackColToEvt: Track fit failed with error code " << error << " track dropped. Number of hits = "<< trkHits.size() << endmsg;  
      
      //delete Track;      
      continue ;
    }
    
    if( track.getNdf() < 0) {       
      debug() << "FullLDCTrackingAlg::AddTrackColToEvt: Track fit returns " << track.getNdf() << " degress of freedom track dropped. Number of hits = "<< trkHits.size() << endmsg;       
      //delete Track;
      continue ;
    }

    edm4hep::TrackState trkStateIP;
    for(int i=0;i<track.trackStates_size();i++){
      trkStateIP = track.getTrackStates(i);
      if(trkStateIP.location ==1/*lcio::TrackState::AtIP*/) break;
    }
    
    if (trkStateIP.location != 1/*lcio::TrackState::AtIP*/) {
      debug() << "FullLDCTrackingAlg::AddTrackColToEvt: Track fit returns " << track.getNdf() << " degress of freedom track dropped. Number of hits = "<< trkHits.size() << endmsg;
      throw EVENT::Exception( std::string("FullLDCTracking_MarlinTrk::AddTrackColToEvt: trkStateIP pointer == NULL ")  ) ;
    }
         
    if (group != NULL) {
      TrackExtendedVec trkVecGrp = group->getTrackExtendedVec();
      int nGrTRK = int(trkVecGrp.size());
      for (int iGr=0;iGr<nGrTRK;++iGr) {
        TrackExtended * subTrack = trkVecGrp[iGr];
        track.addToTracks(subTrack->getTrack());

        // check if it is a tpc looper ...
        if( UTIL::BitSet32( subTrack->getTrack().getType() )[  lcio::ILDDetID::TPC   ] )  {
          
          //const TrackVec segments = subTrack->getTrack().getTracks();
	  std::vector<edm4hep::ConstTrack> segments;
	  std::copy(subTrack->getTrack().tracks_begin(), subTrack->getTrack().tracks_end(), std::back_inserter(segments));
          if ( segments.empty() == false ) {
            
            for (unsigned iSeg=0;iSeg<segments.size();++iSeg) {
              track.addToTracks(segments[iSeg]);
            }

          }

        }
      }
    }
    
    float d0TrkCand = trkCand->getD0();
    float z0TrkCand = trkCand->getZ0();
    //    float phi0TrkCand = trkCand->getPhi();
    // FIXME, fucd
    int nhits_in_vxd = track.getSubDetectorHitNumbers(0);
    int nhits_in_ftd = track.getSubDetectorHitNumbers(1);
    int nhits_in_sit = track.getSubDetectorHitNumbers(2);
    int nhits_in_tpc = track.getSubDetectorHitNumbers(3);
    int nhits_in_set = track.getSubDetectorHitNumbers(4);
    //int nhits_in_vxd = Track->subdetectorHitNumbers()[ 2 * lcio::ILDDetID::VXD - 2 ];
    //int nhits_in_ftd = Track->subdetectorHitNumbers()[ 2 * lcio::ILDDetID::FTD - 2 ];
    //int nhits_in_sit = Track->subdetectorHitNumbers()[ 2 * lcio::ILDDetID::SIT - 2 ];
    //int nhits_in_tpc = Track->subdetectorHitNumbers()[ 2 * lcio::ILDDetID::TPC - 2 ];
    //int nhits_in_set = Track->subdetectorHitNumbers()[ 2 * lcio::ILDDetID::SET - 2 ];
    
    int nHitsSi = nhits_in_vxd + nhits_in_ftd + nhits_in_sit;
    
    debug() << " Hit numbers for Track "<< track.id() << ": "
	    << " vxd hits = " << nhits_in_vxd
	    << " ftd hits = " << nhits_in_ftd
	    << " sit hits = " << nhits_in_sit
	    << " tpc hits = " << nhits_in_tpc
	    << " set hits = " << nhits_in_set
	    << endmsg;
    
    if (nhits_in_vxd > 0) track.setType( track.getType()| (1<<lcio::ILDDetID::VXD) ) ;
    if (nhits_in_ftd > 0) track.setType( track.getType()| (1<<lcio::ILDDetID::FTD) ) ;
    if (nhits_in_sit > 0) track.setType( track.getType()| (1<<lcio::ILDDetID::SIT) ) ;
    if (nhits_in_tpc > 0) track.setType( track.getType()| (1<<lcio::ILDDetID::TPC) ) ;
    if (nhits_in_set > 0) track.setType( track.getType()| (1<<lcio::ILDDetID::SET) ) ;
    
    bool rejectTrack_onTPCHits = (nhits_in_tpc < _cutOnTPCHits) && (nHitsSi<=0);
    
    bool rejectTrackonSiliconHits = ( (nhits_in_tpc<=0) && (nHitsSi<_cutOnSiHits) );
    bool rejectTrackonImpactParameters =  ( fabs(d0TrkCand) > _d0TrkCut ) || ( fabs(z0TrkCand) > _z0TrkCut );
    
    if ( rejectTrack_onTPCHits || rejectTrackonSiliconHits  || rejectTrackonImpactParameters) {

      debug() << " Track " << trkCand
	      << " rejected : rejectTrack_onTPCHits = " << rejectTrack_onTPCHits
	      << " rejectTrackonSiliconHits " << rejectTrackonSiliconHits
	      << " rejectTrackonImpactParameters " << rejectTrackonImpactParameters
	      << endmsg;
      //delete Track;
    }
    else {
      float omega = trkStateIP.omega;
      float tanLambda = trkStateIP.tanLambda;
      float phi0 = trkStateIP.phi;
      float d0 = trkStateIP.D0;
      float z0 = trkStateIP.Z0;
      
      HelixClass helix;
      helix.Initialize_Canonical(phi0,d0,z0,omega,tanLambda,_bField);
      
      float trkPx = helix.getMomentum()[0];
      float trkPy = helix.getMomentum()[1];
      float trkPz = helix.getMomentum()[2];
      float trkP = sqrt(trkPx*trkPx+trkPy*trkPy+trkPz*trkPz);
      
      eTot += trkP;
      pxTot += trkPx;
      pyTot += trkPy;
      pzTot += trkPz;   
      nTotTracks++;
   
      colTRK->push_back(track);
      debug() << " Add Track to final Collection: ID = " << track.id() << " for trkCand "<< trkCand << endmsg;
    }
  }
  
  // streamlog_out(DEBUG5) << endmsg;
  debug() << "Number of accepted " << _OutputTrackColHdl.fullKey() << " = " << nTotTracks << endmsg;
  debug() << "Total 4-momentum of " << _OutputTrackColHdl.fullKey() << " : E = " << eTot
	  << " Px = " << pxTot
	  << " Py = " << pyTot
	  << " Pz = " << pzTot << endmsg;

}

void FullLDCTrackingAlg::prepareVectors() {
  
  _allTPCHits.clear();
  _allVTXHits.clear();
  _allFTDHits.clear();
  _allSITHits.clear();
  _allSETHits.clear();
  _allETDHits.clear();
  _allTPCTracks.clear();
  _allSiTracks.clear();
  _allCombinedTracks.clear();
  _allNonCombinedTPCTracks.clear();
  _allNonCombinedSiTracks.clear();
  _trkImplVec.clear();
  _candidateCombinedTracks.clear();
  
  
  std::map <edm4hep::ConstTrackerHit,TrackerHitExtended*> mapTrackerHits;
  
  // Reading TPC hits
  const edm4hep::TrackerHitCollection* hitTPCCol = nullptr;
  try {
    hitTPCCol = _TPCTrackerHitColHdl.get();
  }
  catch ( GaudiException &e ) {
    debug() << "Collection " << _TPCTrackerHitColHdl.fullKey() << " is unavailable in event " << _nEvt << endmsg;
  }
  if(hitTPCCol){
    int nelem = hitTPCCol->size();
    debug() << "Number of TPC hits = " << nelem << endmsg;
    for (edm4hep::ConstTrackerHit hit : *hitTPCCol) {
      TrackerHitExtended * hitExt = new TrackerHitExtended(hit);
      
      // Covariance Matrix in LCIO is defined in XYZ convert to R-Phi-Z
      // For no error in r
      
      double tpcRPhiRes = sqrt(hit.getCovMatrix()[0] + hit.getCovMatrix()[2]);
      double tpcZRes = sqrt(hit.getCovMatrix()[5]);
      
      hitExt->setResolutionRPhi(float(tpcRPhiRes));
      hitExt->setResolutionZ(float(tpcZRes));
      
      // type and det are no longer used, set to INT_MAX to try and catch any missuse
      hitExt->setType(int(INT_MAX));
      hitExt->setDet(int(INT_MAX));
      _allTPCHits.push_back( hitExt );
      mapTrackerHits[hit] = hitExt;

      double pos[3];
      
      for (int i=0; i<3; ++i) {
        pos[i] = hit.getPosition()[i];
      }
      
      unsigned int layer = static_cast<unsigned int>(getLayerID(hit));
      
      debug() << " TPC Hit added : @ " << pos[0] << " " << pos[1] << " " << pos[2]  << " drphi " << tpcRPhiRes << " dz " << tpcZRes << "  layer = " << layer << endmsg;
    
    }
  }
  
  // Reading in FTD Pixel Hits Collection
  //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  const edm4hep::TrackerHitCollection* hitFTDPixelCol = nullptr;
  try {
    hitFTDPixelCol = _FTDPixelTrackerHitColHdl.get();
  }
  catch ( GaudiException &e ) {
    debug() << "Collection " << _FTDPixelTrackerHitColHdl.fullKey() << " is unavailable in event " << _nEvt << endmsg;
  }

  if(hitFTDPixelCol){
    int nelem = hitFTDPixelCol->size();
    debug() << "Number of FTD Pixel Hits = " << nelem << endmsg;
    for(edm4hep::ConstTrackerHit hit : *hitFTDPixelCol){
      if ( UTIL::BitSet32( hit.getType() )[ UTIL::ILDTrkHitTypeBit::ONE_DIMENSIONAL ] ) continue;

      TrackerHitExtended * hitExt = new TrackerHitExtended( hit );
      
      // double point_res_rphi = sqrt( hit->getdU()*hit->getdU() + hit->getdV()*hit->getdV() );
      // FIXME to calculate the correct resolution;
      double point_res_rphi = sqrt( hit.getCovMatrix()[2]*hit.getCovMatrix()[2] + hit.getCovMatrix()[5]*hit.getCovMatrix()[5] );
      hitExt->setResolutionRPhi( point_res_rphi );
      
      // SJA:FIXME why is this needed? 
      hitExt->setResolutionZ(0.1);
      
      // type and det are no longer used, set to INT_MAX to try and catch any missuse
      hitExt->setType(int(INT_MAX));            
      hitExt->setDet(int(INT_MAX));
      
      _allFTDHits.push_back( hitExt );
      mapTrackerHits[hit] = hitExt;
      
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
        
      }
      
      if (layer >= _nLayersFTD) {
        fatal() << "FullLDCTrackingAlg => fatal error in FTD : layer is outside allowed range : " << layer << " number of layers = " << _nLayersFTD << endmsg;
        exit(1);
      }
      
      double pos[3];
      
      for (int i=0; i<3; ++i) {
        pos[i] = hit.getPosition()[i];
      }      
      
      debug() << " FTD Pixel Hit added : @ " << pos[0] << " " << pos[1] << " " << pos[2]  << " drphi " << hitExt->getResolutionRPhi() << " dz " << hitExt->getResolutionZ() << "  layer = " << layer << endmsg;
      
    }
  }
  
  // Reading in FTD SpacePoint Collection
  //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  const edm4hep::TrackerHitCollection* hitFTDSpacePointCol = nullptr;
  try {
    hitFTDSpacePointCol = _FTDSpacePointColHdl.get();
  }
  catch ( GaudiException &e ) {
    debug() << "Collection " << _FTDSpacePointColHdl.fullKey() << " is unavailable in event " << _nEvt << endmsg;
  }

  const edm4hep::TrackerHitCollection* rawFTDCol = nullptr;
  if(hitFTDSpacePointCol){
    try{
      rawFTDCol = _inFTDRawColHdl.get();
    }
    catch ( GaudiException &e ) {
      fatal() << "Collection " << _inFTDRawColHdl.fullKey() << " is unavailable in event " << _nEvt << endmsg;
    }
  }

  if(hitFTDSpacePointCol&&rawFTDCol){
    Navigation::Instance()->AddTrackerHitCollection(rawFTDCol);
    int nelem = hitFTDSpacePointCol->size();
    debug() << "Number of FTD SpacePoints = " << nelem << endmsg;
    for(edm4hep::ConstTrackerHit hit : *hitFTDSpacePointCol){
      TrackerHitExtended * hitExt = new TrackerHitExtended( hit );
      
      // SJA:FIXME: fudge for now by a factor of two and ignore covariance
      double point_res_rphi = 2 * sqrt( hit.getCovMatrix()[0] + hit.getCovMatrix()[2] );
      
      hitExt->setResolutionRPhi( point_res_rphi );
      
      // SJA:FIXME why is this needed? 
      hitExt->setResolutionZ(0.1);
      
      // type is now only used in one place where it is set to 0 to reject hits from a fit, set to INT_MAX to try and catch any missuse
      hitExt->setType(int(INT_MAX));
      // det is no longer used set to INT_MAX to try and catch any missuse
      hitExt->setDet(int(INT_MAX));
      
      _allFTDHits.push_back( hitExt );
      mapTrackerHits[hit] = hitExt;
      
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
        
      }
      
      if (layer >= _nLayersFTD) {
        fatal() << "FullLDCTrackingAlg => fatal error in FTD : layer is outside allowed range : " << layer << " number of layers = " << _nLayersFTD <<  endmsg;
        exit(1);
      }
      
      double pos[3];
      
      for (int i=0; i<3; ++i) {
        pos[i] = hit.getPosition()[i];
      }      
      
      debug() << " FTD SpacePoint Hit added : @ " << pos[0] << " " << pos[1] << " " << pos[2]  << " drphi " << hitExt->getResolutionRPhi() << " dz " << hitExt->getResolutionZ() << "  layer = " << layer << endmsg;
      
      
    }
  }
  
  const edm4hep::TrackerHitCollection* hitSITCol = nullptr;
  try {
    hitSITCol = _SITTrackerHitColHdl.get();
  }
  catch ( GaudiException &e ) {
    debug() << "Collection " << _SITTrackerHitColHdl.fullKey() << " is unavailable in event " << _nEvt << endmsg;
  }
  
  const edm4hep::TrackerHitCollection* rawSITCol = nullptr;
  if(hitSITCol){
    try{
      rawSITCol = _inSITRawColHdl.get();
    }
    catch ( GaudiException &e ) {
      warning() << "Collection " << _inSITRawColHdl.fullKey() << " is unavailable in event " << _nEvt << ", if SIT is Space Point, it needed " << endmsg;
    }
  }
  
  if(hitSITCol){
    if(rawSITCol) Navigation::Instance()->AddTrackerHitCollection(rawSITCol);
    
    int nelem = hitSITCol->size();

    debug() << "Number of SIT hits = " << nelem << endmsg;
    
    double drphi(NAN);
    double dz(NAN);
    
    for(edm4hep::ConstTrackerHit trkhit : *hitSITCol){
      
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
      
      int layer = getLayerID(trkhit);
      
      if (layer < 0 || (unsigned)layer >= _nLayersSIT) {
        fatal() << "FullLDCTrackingAlg => fatal error in SIT : layer is outside allowed range : " << layer << endmsg;
        exit(1);
      }
      
      // first check that we have not been given 1D hits by mistake, as they won't work here
      if ( UTIL::BitSet32( trkhit.getType() )[ UTIL::ILDTrkHitTypeBit::ONE_DIMENSIONAL ] ) {
        fatal() << "FullLDCTrackingAlg: SIT Hit cannot be of type UTIL::ILDTrkHitTypeBit::ONE_DIMENSIONAL COMPOSITE SPACEPOINTS must be use instead. \n\n  exit(1) called from file " << __FILE__ << " and line " << __LINE__ << endmsg;
        exit(1);
      } 
      // most likely case: COMPOSITE_SPACEPOINT hits formed from stereo strip hits
      else if ( UTIL::BitSet32( trkhit.getType() )[ UTIL::ILDTrkHitTypeBit::COMPOSITE_SPACEPOINT ] ) {
        
        // SJA:FIXME: fudge for now by a factor of two and ignore covariance
        drphi =  2 * sqrt(trkhit.getCovMatrix()[0] + trkhit.getCovMatrix()[2]);         
        dz    =      sqrt(trkhit.getCovMatrix()[5]);         
        
      } 
      // or a PIXEL based SIT, using 2D TrackerHitPlane like the VXD above
      else if ( UTIL::BitSet32( trkhit.getType() )[3])  {
        // FIXME Should make it correct
	// first we need to check if the measurement vectors are aligned with the global coordinates 
	gear::Vector3D U(1.0,trkhit.getCovMatrix()[1],trkhit.getCovMatrix()[0],gear::Vector3D::spherical);
	gear::Vector3D V(1.0,trkhit.getCovMatrix()[4],trkhit.getCovMatrix()[3],gear::Vector3D::spherical);
	gear::Vector3D Z(0.0,0.0,1.0);
        
        const float eps = 1.0e-07;
        // V must be the global z axis 
        if( fabs(1.0 - V.dot(Z)) > eps ) {
          error() << "FullLDCTrackingAlg: PIXEL SIT Hit measurment vectors V is not equal to the global Z axis. \n\n  exit(1) called from file " << __FILE__ << " and line " << __LINE__ << endmsg;
          exit(1);
        }
        
        // U must be normal to the global z axis
        if( fabs(U.dot(Z)) > eps ) {
          error() << "FullLDCTrackingAlg: PIXEL SIT Hit measurment vectors U is not in the global X-Y plane. \n\n exit(1) called from file " << __FILE__ << " and line " << __LINE__ << endmsg;
          exit(1);
        }
        
        // FIXME should make it correct
        // drphi = trkhit_P->getdU();
        // dz    = trkhit_P->getdV();                                                 
        drphi = trkhit.getCovMatrix()[2];
	dz    = trkhit.getCovMatrix()[5];
      } 
      // or a simple cylindrical design, as used in the LOI
      /* FIXME, fucd
      else if ( true ) {
        trkhit_C = hitCollection->at( ielem );
        // FIXME
        // drphi = trkhit_C->getdRPhi();
        // dz    = trkhit_C->getdZ();
        drphi = 1.0;
        dz = 1.0; 
      } 
      */
      // this would be very unlikely, but who knows ... just an ordinary TrackerHit, which is not a COMPOSITE_SPACEPOINT
      else {
        
        // SJA:FIXME: fudge for now by a factor of two and ignore covariance
        drphi =  2 * sqrt(trkhit.getCovMatrix()[0] + trkhit.getCovMatrix()[2]);         
        dz =     sqrt(trkhit.getCovMatrix()[5]);             
        
      }
      
      // now that the hit type has been established carry on and create a 
      
      TrackerHitExtended * hitExt = new TrackerHitExtended( trkhit );
      
      // SJA:FIXME: just use planar res for now
      hitExt->setResolutionRPhi(drphi);
      hitExt->setResolutionZ(dz);
      
      // set type is now only used in one place where it is set to 0 to reject hits from a fit, set to INT_MAX to try and catch any missuse
      hitExt->setType(int(INT_MAX));
      // det is no longer used set to INT_MAX to try and catch any missuse
      hitExt->setDet(int(INT_MAX));
      
      _allSITHits.push_back( hitExt );
      mapTrackerHits[trkhit] = hitExt;
      
      
      double pos[3];
      
      for (int i=0; i<3; ++i) {
        pos[i] = trkhit.getPosition()[i];
      }
      
      debug() << " SIT Hit " <<  trkhit.id() << " added : @ " << pos[0] << " " << pos[1] << " " << pos[2] << " drphi " << hitExt->getResolutionRPhi() << " dz " << hitExt->getResolutionZ() << "  layer = " << layer << endmsg;
    }
  }
  
  const edm4hep::TrackerHitCollection* hitSETCol = nullptr;
  try {
    hitSETCol = _SETTrackerHitColHdl.get();
  }
  catch ( GaudiException &e ) {
    debug() << "Collection " << _SETTrackerHitColHdl.fullKey() << " is unavailable in event " << _nEvt << endmsg;
  }

  const edm4hep::TrackerHitCollection* rawSETCol = nullptr;
  if(hitSETCol){
    try{
      rawSETCol = _inSETRawColHdl.get();
    }
    catch ( GaudiException &e ) {
      warning() << "Collection " << _inSETRawColHdl.fullKey() << " is unavailable in event " << _nEvt << ", if SIT is Space Point, it needed " << endmsg;
    }
  }

  if(hitSETCol){
    if(rawSETCol) Navigation::Instance()->AddTrackerHitCollection(rawSETCol);

    int nelem = hitSETCol->size();

    debug() << "Number of SET hits = " << nelem << endmsg;

    double drphi(NAN);
    double dz(NAN);

    for(edm4hep::ConstTrackerHit trkhit : *hitSETCol){
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
      int layer = getLayerID(trkhit);
      
      if (layer < 0 || (unsigned)layer >= _nLayersSET) {
        fatal() << "FullLDCTrackingAlg => fatal error in SET : layer is outside allowed range : " << layer << endmsg;
        exit(1);
      }
      
      // first check that we have not been given 1D hits by mistake, as they won't work here
      if ( UTIL::BitSet32( trkhit.getType() )[ UTIL::ILDTrkHitTypeBit::ONE_DIMENSIONAL ] ) {
	fatal() << "SiliconTrackingAlg => fatal error in SIT : layer is outside allowed range : " << layer << endmsg;
        exit(1);
      }
      // most likely case: COMPOSITE_SPACEPOINT hits formed from stereo strip hits
      else if ( UTIL::BitSet32( trkhit.getType() )[ UTIL::ILDTrkHitTypeBit::COMPOSITE_SPACEPOINT ] ) {
	// SJA:FIXME: fudge for now by a factor of two and ignore covariance
        drphi =  2 * sqrt(trkhit.getCovMatrix()[0] + trkhit.getCovMatrix()[2]);
        dz    =      sqrt(trkhit.getCovMatrix()[5]);
      }
      // or a PIXEL based SET, using 2D TrackerHitPlane like the VXD above
      else if ( UTIL::BitSet32( trkhit.getType() )[3] )  {
        // FIXME should consider this carefully
        
        // first we need to check if the measurement vectors are aligned with the global coordinates
        // gear::Vector3D U(1.0,trkhit_P->getU()[1],trkhit_P->getU()[0],gear::Vector3D::spherical);
        // gear::Vector3D V(1.0,trkhit_P->getV()[1],trkhit_P->getV()[0],gear::Vector3D::spherical);
        // FIXME Should calculate it correctly
	gear::Vector3D U(1.0,trkhit.getCovMatrix()[1],trkhit.getCovMatrix()[0],gear::Vector3D::spherical);
	gear::Vector3D V(1.0,trkhit.getCovMatrix()[4],trkhit.getCovMatrix()[3],gear::Vector3D::spherical);
	gear::Vector3D Z(0.0,0.0,1.0);
        
        const float eps = 1.0e-07;
        // V must be the global z axis
        if( fabs(1.0 - V.dot(Z)) > eps ) {
          fatal() << "FullLDCTrackingAlg: PIXEL SET Hit measurment vectors V is not equal to the global Z axis. \n\n  exit(1) called from file " << __FILE__ << " and line " << __LINE__ << endmsg;
          exit(1);
        }
        
        // U must be normal to the global z axis
        if( fabs(U.dot(Z)) > eps ) {
          fatal() << "FullLDCTrackingAlg: PIXEL SET Hit measurment vectors U is not in the global X-Y plane. \n\n exit(1) called from file " << __FILE__ << " and line " << __LINE__ << endmsg;
          exit(1);
        }
	
        // FIXME should use the correct 
        // drphi = trkhit_P->getdU();
        // dz    = trkhit_P->getdV();
        drphi = trkhit.getCovMatrix()[2];
	dz    = trkhit.getCovMatrix()[5];
      }
      // or a simple cylindrical design, as used in the LOI
      /*FIXME, fucd
      else if ( true ) {
        trkhit_C = hitCollection->at( ielem );
        // FIXME
        // drphi = trkhit_C->getdRPhi();
        // dz    = trkhit_C->getdZ();
        drphi = 1.0;
        dz = 1.0;
      }
      */
      // this would be very unlikely, but who knows ... just an ordinary TrackerHit, which is not a COMPOSITE_SPACEPOINT
      else {
        
        // SJA:FIXME: fudge for now by a factor of two and ignore covariance
        drphi =  2 * sqrt(trkhit.getCovMatrix()[0] + trkhit.getCovMatrix()[2]);
        dz =     sqrt(trkhit.getCovMatrix()[5]);
      }
      
      // now that the hit type has been established carry on and create a
      
      TrackerHitExtended * hitExt = new TrackerHitExtended( trkhit );
      
      // SJA:FIXME: just use planar res for now
      hitExt->setResolutionRPhi(drphi);
      hitExt->setResolutionZ(dz);
      
      // set type is now only used in one place where it is set to 0 to reject hits from a fit, set to INT_MAX to try and catch any missuse
      hitExt->setType(int(INT_MAX));
      // det is no longer used set to INT_MAX to try and catch any missuse
      hitExt->setDet(int(INT_MAX));
      
      _allSETHits.push_back( hitExt );
      mapTrackerHits[trkhit] = hitExt;
      
      
      double pos[3];
      
      for (int i=0; i<3; ++i) {
        pos[i] = trkhit.getPosition()[i];
      }
      
      debug() << " SET Hit " <<  trkhit.id() << " added : @ " << pos[0] << " " << pos[1] << " " << pos[2] << " drphi " << hitExt->getResolutionRPhi() << " dz " << hitExt->getResolutionZ() << "  layer = " << layer << endmsg;
    }
  }
    
  // Reading VTX Hits
  const edm4hep::TrackerHitCollection* hitVTXCol = nullptr;
  try {
    hitVTXCol = _VTXTrackerHitColHdl.get();
  }
  catch ( GaudiException &e ) {
    debug() << "Collection " << _VTXTrackerHitColHdl.fullKey() << " is unavailable in event " << _nEvt << endmsg;
  }

  if(hitVTXCol){
    int nelem = hitVTXCol->size();
    debug() << "Number of VTX hits = " << nelem << endmsg;
    for(edm4hep::ConstTrackerHit trkhit : *hitVTXCol){
      // FIXME tracker hit plane type of the hit
      TrackerHitExtended* hitExt = new TrackerHitExtended(trkhit);
      
      // SJA:FIXME: just use planar res for now
      // FIXME should use the correct resolution
      // hitExt->setResolutionRPhi(trkhit->getdU());
      // hitExt->setResolutionZ(trkhit->getdV());
      hitExt->setResolutionRPhi(trkhit.getCovMatrix()[2]);
      hitExt->setResolutionZ(trkhit.getCovMatrix()[5]);

      // type and det are no longer used, set to INT_MAX to try and catch any missuse
      hitExt->setType(int(INT_MAX));      
      hitExt->setDet(int(INT_MAX));
      _allVTXHits.push_back( hitExt );
      mapTrackerHits[trkhit] = hitExt;
      
      double pos[3];
      
      for (int i=0; i<3; ++i) {
        pos[i] = trkhit.getPosition()[i];
      }
      
      int layer = getLayerID(trkhit);
      
      debug() << " VXD Hit " <<  trkhit.id() << " added : @ " << pos[0] << " " << pos[1] << " " << pos[2] << " drphi " << hitExt->getResolutionRPhi() << " dz " << hitExt->getResolutionZ() << "  layer = " << layer << endmsg;
    }
  }
    
  // Reading TPC Tracks
  const edm4hep::TrackCollection* tpcTrackCol = nullptr;
  try {
    tpcTrackCol = _TPCTrackColHdl.get();
  }
  catch ( GaudiException &e ) {
    debug() << "Collection " << _TPCTrackColHdl.fullKey() << " is unavailable in event " << _nEvt << endmsg;
  }
  if(tpcTrackCol){
    int nelem = tpcTrackCol->size();
    debug() << "Number of TPC Tracks = " << nelem << endmsg;
    debug() << " Trk    ID        p          D0         Z0       Px       Py       Pz   ntpc ndf  Chi2/ndf nlinkedTracks" << endmsg;
    int iTrk = -1;
    for(edm4hep::ConstTrack tpcTrack : *tpcTrackCol){
      iTrk++;
      ConstTrackerHitVec hitVec(tpcTrack.trackerHits_begin(), tpcTrack.trackerHits_end());
      int nHits = int(hitVec.size());
      
      debug() << toString( iTrk, tpcTrack ,  _bField ) << endmsg;
      
      TrackExtended * trackExt = new TrackExtended( tpcTrack );

      trackExt->setOmega(getOmega(tpcTrack));
      trackExt->setTanLambda(getTanLambda(tpcTrack));
      trackExt->setPhi(getPhi(tpcTrack));
      trackExt->setD0(getD0(tpcTrack));
      trackExt->setZ0(getZ0(tpcTrack));
      float cov[15];
      float param[5];
      //      float reso[4];
      param[0] = getOmega(tpcTrack);
      param[1] = getTanLambda(tpcTrack);
      param[2] = getPhi(tpcTrack);
      param[3] = getD0(tpcTrack);
      param[4] = getZ0(tpcTrack);
      
      std::array<float, 15> Cov = getCovMatrix(tpcTrack);
      int NC = int(Cov.size());
      for (int ic=0;ic<NC;ic++) {
        cov[ic] =  Cov[ic];
      }
      
      trackExt->setCovMatrix(cov);
      trackExt->setNDF(tpcTrack.getNdf());
      trackExt->setChi2(tpcTrack.getChi2());            
      for (int iHit=0;iHit<nHits;++iHit) {
	edm4hep::ConstTrackerHit hit = hitVec[iHit];
        TrackerHitExtended * hitExt = mapTrackerHits[hit];
        hitExt->setTrackExtended( trackExt );
        trackExt->addTrackerHitExtended( hitExt );      
      }      
      _allTPCTracks.push_back( trackExt );                
    }      
  }

  // Reading Si Tracks
  const edm4hep::TrackCollection* siTrackCol = nullptr;
  try {
    siTrackCol = _SiTrackColHdl.get();
  }
  catch ( GaudiException &e ) {
    debug() << "Collection " << _SiTrackColHdl.fullKey() << " is unavailable in event " << _nEvt << endmsg;
  }
  if(siTrackCol){
    int nelem = siTrackCol->size();
    debug() << "Number of Si Tracks = " << nelem << endmsg;
    debug() << " Trk    ID        p          D0         Z0       Px       Py       Pz   ntpc ndf  Chi2/ndf nlinkedTracks" << endmsg;
    int iTrk = -1;
    for(edm4hep::ConstTrack siTrack : *siTrackCol){
      iTrk++;
      double prob = ( siTrack.getNdf() > 0 ? gsl_cdf_chisq_Q(  siTrack.getChi2() ,  (double) siTrack.getNdf() )  : 0. ) ;
      if( prob < _minChi2ProbForSiliconTracks ) {
        debug() << "Si Tracks " << siTrack << " id : " << siTrack.id() << " rejected with prob " << prob << " < " << _minChi2ProbForSiliconTracks << endmsg;
        continue;
      }
      
      TrackExtended * trackExt = new TrackExtended( siTrack );
      ConstTrackerHitVec hitVec(siTrack.trackerHits_begin(), siTrack.trackerHits_end());
      int nHits = int(hitVec.size());
      trackExt->setOmega(getOmega(siTrack));
      trackExt->setTanLambda(getTanLambda(siTrack));
      trackExt->setPhi(getPhi(siTrack));
      trackExt->setD0(getD0(siTrack));
      trackExt->setZ0(getZ0(siTrack));
      float cov[15];
      float param[5];
      
      param[0] = getOmega(siTrack);
      param[1] = getTanLambda(siTrack);
      param[2] = getPhi(siTrack);
      param[3] = getD0(siTrack);
      param[4] = getZ0(siTrack);      
            
      std::array<float, 15> Cov = getCovMatrix(siTrack);
      int NC = int(Cov.size());
      for (int ic=0;ic<NC;ic++) {
        cov[ic] =  Cov[ic];
      } 
      trackExt->setCovMatrix(cov);
      trackExt->setNDF(siTrack.getNdf());
      trackExt->setChi2(siTrack.getChi2());      
      char strg[200];
      HelixClass helixSi;
      for (int iHit=0;iHit<nHits;++iHit) {
	edm4hep::ConstTrackerHit hit = hitVec[iHit];
        TrackerHitExtended * hitExt = mapTrackerHits[hit];
        hitExt->setTrackExtended( trackExt );
        
        trackExt->addTrackerHitExtended( hitExt );      
      }
      
      float d0Si = trackExt->getD0();
      float z0Si = trackExt->getZ0();
      float omegaSi = trackExt->getOmega();
      float phi0Si = trackExt->getPhi();
      float tanLSi = trackExt->getTanLambda();
      helixSi.Initialize_Canonical(phi0Si,d0Si,z0Si,omegaSi,tanLSi,_bField);
      float pxSi = helixSi.getMomentum()[0];
      float pySi = helixSi.getMomentum()[1];
      float pzSi = helixSi.getMomentum()[2];
      const float pTot = sqrt(pxSi*pxSi+pySi*pySi+pzSi*pzSi);
      const int ndfSi = trackExt->getNDF();
      float Chi2Si = trackExt->getChi2()/float(trackExt->getNDF());
      sprintf(strg,"%3i   %5i %9.3f  %9.3f  %9.3f  %7.2f  %7.2f  %7.2f %4i %4i %8.3f",iTrk, siTrack.id(), pTot, d0Si,z0Si,pxSi,pySi,pzSi,nHits, ndfSi, Chi2Si);
      debug() << strg << endmsg;
      
      if(nHits>0){
        _allSiTracks.push_back( trackExt );
      }else{
        delete trackExt;
      }
    }
  }
}


void FullLDCTrackingAlg::CleanUp(){
  
  int nNonCombTpc = int(_allNonCombinedTPCTracks.size());  
  for (int i=0;i<nNonCombTpc;++i) {
    TrackExtended * trkExt = _allNonCombinedTPCTracks[i];
    GroupTracks * group = trkExt->getGroupTracks();
    delete group;
  }
  _allNonCombinedTPCTracks.clear();
  
  int nNonCombSi = int(_allNonCombinedSiTracks.size());  
  for (int i=0;i<nNonCombSi;++i) {
    TrackExtended * trkExt = _allNonCombinedSiTracks[i];
    GroupTracks * group = trkExt->getGroupTracks();
    delete group;
  }
  _allNonCombinedSiTracks.clear();
  
  int nSITHits = int(_allSITHits.size());
  for (int i=0;i<nSITHits;++i) {
    TrackerHitExtended * hitExt = _allSITHits[i];
    delete hitExt;
  }
  _allSITHits.clear();
  
  int nSETHits = int(_allSETHits.size());
  for (int i=0;i<nSETHits;++i) {
    TrackerHitExtended * hitExt = _allSETHits[i];
    delete hitExt;
  }
  _allSETHits.clear();
  
  int nTPCHits = int(_allTPCHits.size());
  for (int i=0;i<nTPCHits;++i) {
    TrackerHitExtended * hitExt = _allTPCHits[i];
    delete hitExt;
  }
  _allTPCHits.clear();
  
  int nFTDHits = int(_allFTDHits.size());
  for (int i=0;i<nFTDHits;++i) {
    TrackerHitExtended * hitExt = _allFTDHits[i];
    delete hitExt;
  }
  _allFTDHits.clear();
  
  int nETDHits = int(_allETDHits.size());
  for (int i=0;i<nETDHits;++i) {
    TrackerHitExtended * hitExt = _allETDHits[i];
    delete hitExt;
  }
  _allETDHits.clear();
  
  int nVTXHits = int(_allVTXHits.size());
  for (int i=0;i<nVTXHits;++i) {
    TrackerHitExtended * hitExt = _allVTXHits[i];
    delete hitExt;
  }
  _allVTXHits.clear();
  
  int nSiTrk = int(_allSiTracks.size());
  for (int i=0;i<nSiTrk;++i) {
    TrackExtended * trkExt = _allSiTracks[i];
    delete trkExt;
  }
  _allSiTracks.clear();
  
  int nTPCTrk = int(_allTPCTracks.size());
  for (int i=0;i<nTPCTrk;++i) {
    TrackExtended * trkExt = _allTPCTracks[i];
    delete trkExt;
  }
  _allTPCTracks.clear();
  
  int nCombTrk = int(_allCombinedTracks.size());
  for (int i=0;i<nCombTrk;++i) {
    TrackExtended * trkExt = _allCombinedTracks[i];
    GroupTracks * group = trkExt->getGroupTracks();
    delete group;
    delete trkExt;    
  }
  _allCombinedTracks.clear();
  
  //   int nImplTrk = int(_trkImplVec.size());
  //   for (int i=0;i<nImplTrk;++i) {
  //     TrackExtended * trkImpl = _trkImplVec[i];
  //     delete trkImpl;
  //   }
  _trkImplVec.clear();
  
  //AS: Dont delete the individual entries, some of them are cleared elsewhere, I think
  _candidateCombinedTracks.clear();
}

/*
 
 Compare all TPC tracks with all Silicon tracks
 and compare using CompareTrkII which really only compares delta d0 and delta z0
 then calculate dOmega and the angle between the tracks, which are checked against cuts here.
 
 if CompareTrkII and the cuts on dOmega and angle succecede try and full fit of the hits using CombineTracks
 assuming that the merger was not vetoed VetoMerge
 i.e. if the momentum of either track is less than 2.5 GeV
 or if following a full fit the NDF+10 of the combined tracks is less than the NDF_first + NDF_second
 
 if this is successful the track is added to _allCombinedTracks and the TPC and Si Segements are added to _candidateCombinedTracks
 
 */

void FullLDCTrackingAlg::MergeTPCandSiTracks() {
  
  int nTPCTracks = int(_allTPCTracks.size());
  int nSiTracks  = int(_allSiTracks.size());
  
  debug() << " MergeTPCandSiTracks called nTPC tracks " << nTPCTracks << " - nSiTracks " << nSiTracks << endmsg;
  
  for (int iTPC=0;iTPC<nTPCTracks;++iTPC) {
    TrackExtended * tpcTrackExt = _allTPCTracks[iTPC];
    for (int iSi=0;iSi<nSiTracks;++iSi) {
      TrackExtended * siTrackExt = _allSiTracks[iSi];
      int iComp = 0;
      float angle = 0;

      debug() << " compare tpc trk " << toString( iTPC,  tpcTrackExt->getTrack(), _bField  ) << endmsg;
      debug() << "    to si trk    " << toString( iSi,   siTrackExt->getTrack(),  _bField  ) << endmsg;
      
      float dOmega = CompareTrkII(siTrackExt,tpcTrackExt,_d0CutForMerging,_z0CutForMerging,iComp,angle);
      
      if ( (dOmega<_dOmegaForMerging) && (angle<_angleForMerging) && !VetoMerge(tpcTrackExt,siTrackExt)) {
	
        debug() << " call CombineTracks for tpc trk " << tpcTrackExt << " si trk " << siTrackExt << endmsg;
	
        TrackExtended *combinedTrack = CombineTracks(tpcTrackExt,siTrackExt,_maxAllowedPercentageOfOutliersForTrackCombination, false);
	
        debug() << " combinedTrack returns " << combinedTrack << endmsg;
        
        if (combinedTrack != NULL) {


          _allCombinedTracks.push_back( combinedTrack );
          _candidateCombinedTracks.insert(tpcTrackExt);
          _candidateCombinedTracks.insert(siTrackExt);

	  //          streamlog_out(DEBUG3) << " combinedTrack successfully added to _allCombinedTracks : " << toString( 0,combinedTrack->getTrack(), _bField  )  << endmsg;
          debug() << " *** combinedTrack successfully added to _allCombinedTracks : tpc " << iTPC << " si " << iSi   << endmsg;
	  
          if (_debug >= 3 ) {
            int iopt = 1;
            // PrintOutMerging(tpcTrackExt,siTrackExt,iopt);
          }
        }
	else{
          if (_debug >= 3 ) {
            int iopt = 6;
            // PrintOutMerging(tpcTrackExt,siTrackExt,iopt);
          }
        }
      }
      else {
        if (_debug >= 3 ) {
          int iopt = 6;
          // PrintOutMerging(tpcTrackExt,siTrackExt,iopt);
        }
      }
    }
  }
  
  
}


/*
 
 All Si and TPC tracks which have been merged using MergeTPCandSiTracks are excluded from this search.
 
 The remaining TPC tracks and Silicon tracks will be tested using CompareTrkIII 

 CompareTrkIII does the following 
 
 i)   significance d0    < d0Cut
 ii)  significance z0    < z0Cut
 iii) pdot > 0.999
 
 if above three cuts are passed return:
 
 the significance dAngle (return by reference)
 and
 the significance dOmega

 if CompareTrkIII and the cuts on the significance of dOmega and angle succecede try and full fit of the hits using CombineTracks
 assuming that the merger was not vetoed VetoMerge 
 i.e. if the momentum of either track is less than 2.5 GeV
      and there is no overlapping of the segments
 
 if this is successful the track is added to _allCombinedTracks and the TPC and Si Segements are added to _candidateCombinedTracks
 
 */

void FullLDCTrackingAlg::MergeTPCandSiTracksII() {
  
  int nTPCTracks = int(_allTPCTracks.size());
  int nSiTracks  = int(_allSiTracks.size());
  
  debug() << " MergeTPCandSiTracksII called nTPC tracks " << nTPCTracks << " - nSiTracks " << nSiTracks << endmsg;

  for (int iTPC=0;iTPC<nTPCTracks;++iTPC) {
    
    // check if the tpc track has already been merged with CompareTrkII
    TrackExtended * tpcTrackExt = _allTPCTracks[iTPC];
    if(_candidateCombinedTracks.find(tpcTrackExt) != _candidateCombinedTracks.end() )continue;
    
    for (int iSi=0;iSi<nSiTracks;++iSi) {
      
      // check if the tpc track has already been merged with CompareTrkII
      TrackExtended * siTrackExt = _allSiTracks[iSi];
      if(_candidateCombinedTracks.find(siTrackExt)!= _candidateCombinedTracks.end() )continue;

      int iComp = 0;
      float angleSignificance = 0;

      float significance = CompareTrkIII(siTrackExt,tpcTrackExt,_d0CutForMerging,_z0CutForMerging,iComp,angleSignificance);

      debug() << " MergeTPCandSiTracksII - tpctrk " << iTPC << " - " << iSi <<  " - significance " << significance
	      << " angleSignificance " << angleSignificance << endmsg;
      
      if ( (significance<10) && (angleSignificance<5) && !VetoMerge(tpcTrackExt,siTrackExt) ) {

        TrackExtended * combinedTrack = CombineTracks(tpcTrackExt,siTrackExt,_maxAllowedPercentageOfOutliersForTrackCombination, false);
        
        debug() << " combinedTrack returns " << combinedTrack << endmsg;
        
	if (combinedTrack != NULL) {
          
          _allCombinedTracks.push_back( combinedTrack );
	  debug() << " *** combinedTrack successfully added to _allCombinedTracks : tpc " << iTPC << " si " << iSi   << endmsg;
	  
	  if (_debug >= 3 ) {
            int iopt = 1;
            // PrintOutMerging(tpcTrackExt,siTrackExt,iopt);
          }
        }
	else{
          if (_debug >= 3 ) {
            int iopt = 6;
            // PrintOutMerging(tpcTrackExt,siTrackExt,iopt);
          }
        }
      }
      else {
        if (_debug >= 3 ) {
          int iopt = 6;
          // PrintOutMerging(tpcTrackExt,siTrackExt,iopt);
        }
      }
    }
  }
}

// if testCombinationOnly is true then hits will not be assigned to the tracks 
TrackExtended * FullLDCTrackingAlg::CombineTracks(TrackExtended * tpcTrack, TrackExtended * siTrack, float maxAllowedOutliers, bool testCombinationOnly) {
  
  TrackExtended * OutputTrack = NULL;
  
  TrackerHitExtendedVec siHitVec = siTrack->getTrackerHitExtendedVec();
  TrackerHitExtendedVec tpcHitVec = tpcTrack->getTrackerHitExtendedVec();
  
  int nSiHits = int(siHitVec.size());
  int nTPCHits = int(tpcHitVec.size());
  int nHits = nTPCHits + nSiHits;
  
  //std::cout << "FullLDCTrackingAlg::CombineTracks nSiHits = " << nSiHits << endmsg;
  //std::cout << "FullLDCTrackingAlg::CombineTracks nTPCHits = " << nTPCHits << endmsg;
  
  ConstTrackerHitVec trkHits;
  trkHits.reserve(nHits);
  
  for (int ih=0;ih<nSiHits;++ih) {
    edm4hep::ConstTrackerHit trkHit = siHitVec[ih]->getTrackerHit();
    if(trkHit.isAvailable()) { 
      trkHits.push_back(trkHit);   
    }
    else{
      throw EVENT::Exception( std::string("FullLDCTrackingAlg::CombineTracks: TrackerHit pointer == NULL ")  ) ;
    }
  }  
  
  for (int ih=0;ih<nTPCHits;++ih) {
    edm4hep::ConstTrackerHit trkHit = tpcHitVec[ih]->getTrackerHit();
    if(trkHit.isAvailable()) { 
      trkHits.push_back(trkHit);   
    }
    else{
      throw EVENT::Exception( std::string("FullLDCTrackingAlg::CombineTracks: TrackerHit pointer == NULL ")  ) ;
    }
  }      
  
  double chi2_D;
  int ndf;
  
  if( trkHits.size() < 3 ) { 
    
    return 0 ;
    
  }
  
  debug() << "FullLDCTrackingAlg::CombineTracks: Sorting Hits " << trkHits.size() << endmsg;
  
  std::vector< std::pair<float, edm4hep::ConstTrackerHit> > r2_values;
  r2_values.reserve(trkHits.size());
  
  for (ConstTrackerHitVec::iterator it=trkHits.begin(); it!=trkHits.end(); ++it) {
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
  
  debug() << "FullLDCTrackingAlg::CombineTracks: Start Fitting: AddHits: number of hits to fit " << trkHits.size() << endmsg;
  
  std::auto_ptr<MarlinTrk::IMarlinTrack> marlin_trk_autop(_trksystem->createTrack());
  MarlinTrk::IMarlinTrack& marlin_trk = *marlin_trk_autop.get();
  
  edm4hep::TrackState pre_fit ;

  int errorCode = IMarlinTrack::success;

  try{
    pre_fit = getTrackStateAt(tpcTrack->getTrack(), 3/*lcio::TrackState::AtLastHit*/);
  }
  catch(std::runtime_error& e){
    error() << e.what() << " should be checked (TPC track)" << endmsg;
    try{
      pre_fit = getTrackStateAt(siTrack->getTrack(), 3/*lcio::TrackState::AtLastHit*/);
    }
    catch(std::runtime_error& e){
      error() << e.what() << " should be checked (Si track)" << endmsg;
    }
  }
  
  // setup initial dummy covariance matrix
  std::array<float,15> covMatrix;
  
  for (unsigned icov = 0; icov<covMatrix.size(); ++icov) {
    covMatrix[icov] = 0;
  }
  
  covMatrix[0]  = ( _initialTrackError_d0    ); //sigma_d0^2
  covMatrix[2]  = ( _initialTrackError_phi0  ); //sigma_phi0^2
  covMatrix[5]  = ( _initialTrackError_omega ); //sigma_omega^2
  covMatrix[9]  = ( _initialTrackError_z0    ); //sigma_z0^2
  covMatrix[14] = ( _initialTrackError_tanL  ); //sigma_tanl^2
  
  pre_fit.covMatrix = covMatrix;
  
  errorCode = MarlinTrk::createFit( trkHits, &marlin_trk, &pre_fit, _bField, IMarlinTrack::backward , _maxChi2PerHit );
  
  if ( errorCode != IMarlinTrack::success ) {
    debug() << "FullLDCTrackingAlg::CombineTracks: creation of fit fails with error " << errorCode << endmsg;
    return 0;
  }
  
  edm4hep::Vector3d point(0.,0.,0.); // nominal IP
  
  edm4hep::TrackState trkState ;
  errorCode = marlin_trk.propagate(point, trkState, chi2_D, ndf ) ;
  
  if ( errorCode != IMarlinTrack::success ) {
    debug() << "FullLDCTrackingAlg::CombineTracks: propagate to IP fails with error " << errorCode << endmsg;    
    return 0;
  }
  
  if ( ndf < 0  ) {
    debug() << "FullLDCTrackingAlg::CombineTracks: Fit failed NDF is less that zero  " << ndf << endmsg;
    return 0;
  }
  
  float chi2Fit = chi2_D/float(ndf);
  
  if ( chi2Fit > _chi2FitCut ) {
    debug() << "FullLDCTrackingAlg::CombineTracks: track fail Chi2 cut of " << _chi2FitCut << " chi2 of track = " <<  chi2Fit << endmsg;
    return 0;
  }
    
  debug() << "FullLDCTrackingAlg::CombineTracks: Check for outliers " << endmsg;
  
  std::vector<std::pair<edm4hep::ConstTrackerHit, double> > outliers ;
  marlin_trk.getOutliers(outliers);
  
  float outlier_pct = outliers.size()/float(trkHits.size()) ;
  
  if ( outlier_pct > maxAllowedOutliers) {
    debug() << "FullLDCTrackingAlg::CombineTracks: percentage of outliers " << outlier_pct << " is greater than cut maximum: " << maxAllowedOutliers << endmsg;
    return 0;
  }

  // sort the hits into outliers from TPC and Silicon as we will reject the combination if more that 2 Si hits get rejected ...
  
  std::vector<TrackerHitExtended*> siHitInFit;
  std::vector<TrackerHitExtended*> siOutliers;
  std::vector<TrackerHitExtended*> tpcHitInFit;
  std::vector<TrackerHitExtended*> tpcOutliers;
  
  for (int i=0;i<nSiHits;++i) {
    
    bool hit_is_outlier = false;
    
    // we need to make sure that in the case of a composite hit we reject this as well
    ConstTrackerHitVec hits;
    
    // all hits, both the 2D tracker hit, as well as any raw hits which belong to it
    ConstTrackerHit hit = siHitVec[i]->getTrackerHit();
    hits.push_back(hit);
    
    // add the raw hits ...
    int nRawHits = hit.rawHits_size();
    if ( nRawHits>0 ) {
      for (unsigned ihit=0; ihit < nRawHits; ++ihit) {
	try{
	  int type = hit.getType();
	  if(UTIL::BitSet32(type)[UTIL::ILDTrkHitTypeBit::COMPOSITE_SPACEPOINT]){
	    edm4hep::ConstTrackerHit rawHit = Navigation::Instance()->GetTrackerHit(hit.getRawHits(ihit));
	    hits.push_back(rawHit);
	  }
	  else debug() << "not space point, id=" << hit.id() << endmsg; 
	}
	catch(std::runtime_error& e){
	  debug() << e.what() << " raw hit id = " << hit.getRawHits(ihit) << " for TrackerHit id = " << hit.id() << endmsg;
	}
      }
    }
    
    // now double loop over the outliers and the hits assosiated with this TrackerHitExtended and compare
    for ( unsigned ohit = 0; ohit < outliers.size(); ++ohit) {
      for (unsigned ihit = 0; ihit < hits.size(); ++ihit) {
        
        // compare outlier pointer to TrackerHit pointer
        if( outliers[ohit].first == hits[ihit] ){
          // silicon outlier found so add the TrackerHitExtended to the list of outliers
          hit_is_outlier = true;
          siOutliers.push_back(siHitVec[i]);
          break;
        }
        if (hit_is_outlier) {
          break;
        }
      }
    }
    
    if( hit_is_outlier == false ){
      // add the TrackerHitExtended to the list of silicon hits used in the fit
      siHitInFit.push_back(siHitVec[i]);
    }
    
  }
  
  // more simple as TPC hits are never composite
  for (int i=0;i<nTPCHits;++i) {
    
    bool hit_is_outlier = false;
    
    for ( unsigned ihit = 0; ihit < outliers.size(); ++ihit) {

      // compare outlier pointer to TrackerHit pointer
      if( outliers[ihit].first == tpcHitVec[i]->getTrackerHit() ){
        // tpc outlier found so add the TrackerHitExtended to the list of outliers
        hit_is_outlier = true;
        tpcOutliers.push_back(tpcHitVec[i]);
        break;
      }
    }
    
    if( hit_is_outlier == false ){
      // add the TrackerHitExtended to the list of tpc hits used in the fit
      tpcHitInFit.push_back(tpcHitVec[i]);
    }
  }
  
  debug() << "FullLDCTrackingAlg::CombineTracks: Check for Silicon Hit rejections ... " << endmsg;
  
  if ( (int)siOutliers.size() > _maxAllowedSiHitRejectionsForTrackCombination ) {
    debug() << "FullLDCTrackingAlg::CombineTracks: Fit rejects " << siOutliers.size() << " silicon hits : max allowed rejections = " << _maxAllowedSiHitRejectionsForTrackCombination << " : Combination rejected " << endmsg;
    return 0;
  }
  
  float omega = trkState.omega;
  float tanlambda = trkState.tanLambda;
  float phi0 = trkState.phi;
  float d0 = trkState.D0;
  float z0 = trkState.Z0;
  
  OutputTrack = new TrackExtended();

  GroupTracks * group = new GroupTracks();
  OutputTrack->setGroupTracks(group);
  
  group->addTrackExtended(siTrack);
  group->addTrackExtended(tpcTrack);
  
  // note OutputTrack which is of type TrackExtended, only takes fits set for ref point = 0,0,0
  OutputTrack->setOmega(omega);
  OutputTrack->setTanLambda(tanlambda);
  OutputTrack->setPhi(phi0);
  OutputTrack->setZ0(z0);
  OutputTrack->setD0(d0);
  OutputTrack->setChi2(chi2_D);
  OutputTrack->setNDF(ndf);
  
  float cov[15];
  
  for (int i = 0 ; i<15 ; ++i) {
    cov[i] = trkState.covMatrix[i];
  }
  
  OutputTrack->setCovMatrix(cov);
  
  // if this is not just a test of the combination add the hits to the combined track
  if ( testCombinationOnly == false ) {

    for (unsigned ihit=0; ihit<siHitInFit.size(); ++ihit) {
      TrackerHitExtended * hitExt = siHitInFit[ihit];
      OutputTrack->addTrackerHitExtended(hitExt);
      hitExt->setUsedInFit(true);
    }
    
    for (unsigned ihit=0; ihit<siOutliers.size(); ++ihit) {
      TrackerHitExtended * hitExt = siOutliers[ihit];
      OutputTrack->addTrackerHitExtended(hitExt);
      hitExt->setUsedInFit(false);
    }

    for (unsigned ihit=0; ihit<tpcHitInFit.size(); ++ihit) {
      TrackerHitExtended * hitExt = tpcHitInFit[ihit];
      OutputTrack->addTrackerHitExtended(hitExt);
      hitExt->setUsedInFit(true);
    }
    
    for (unsigned ihit=0; ihit<tpcOutliers.size(); ++ihit) {
      TrackerHitExtended * hitExt = tpcOutliers[ihit];
      OutputTrack->addTrackerHitExtended(hitExt);
      hitExt->setUsedInFit(false);
    }
  }
    
  debug() << "FullLDCTrackingAlg::CombineTracks: merged track created  " << OutputTrack << " with " << OutputTrack->getTrackerHitExtendedVec().size() << " hits, nhits tpc " << nTPCHits << " nSiHits " << nSiHits << ", testCombinationOnly = " << testCombinationOnly << endmsg;
  
  return OutputTrack;
  
}


void FullLDCTrackingAlg::SortingTrackHitPairs(TrackHitPairVec & trackHitPairVec) {
  
  int sizeOfVector = int(trackHitPairVec.size());
  TrackHitPair *one,*two,*Temp;
  
  for (int i = 0 ; i < sizeOfVector-1; i++)
    for (int j = 0; j < sizeOfVector-i-1; j++) {
      one = trackHitPairVec[j];
      two = trackHitPairVec[j+1];
      float oneQ = one->getDistance();
      float twoQ = two->getDistance();
      if( oneQ > twoQ ) {
        Temp = trackHitPairVec[j];
        trackHitPairVec[j] = trackHitPairVec[j+1];
        trackHitPairVec[j+1] = Temp;
      }
    }  
  
  
}

/*
 
 Sorts all tracks in the vector by Chi2/NDF
 
 */

void FullLDCTrackingAlg::Sorting(TrackExtendedVec & trackVec) {
  
  int sizeOfVector = int(trackVec.size());
  TrackExtended *one,*two,*Temp;
  
  for (int i = 0 ; i < sizeOfVector-1; i++)
    for (int j = 0; j < sizeOfVector-i-1; j++) {
      one = trackVec[j];
      two = trackVec[j+1];
      float oneQ = one->getChi2()/float(one->getNDF());
      float twoQ = two->getChi2()/float(two->getNDF());
      if( oneQ > twoQ ) {
        Temp = trackVec[j];
        trackVec[j] = trackVec[j+1];
        trackVec[j+1] = Temp;
      }
    }  
}

/*
 
 this function is used to select only the combined tracks which have been formed from excatly 2 sub tracks which themselves are not formed from any other tracks.
 
 */

void FullLDCTrackingAlg::SelectCombinedTracks() {
  
  int nCombTrk = int(_allCombinedTracks.size());
  
  debug() << " **SelectCombinedTracks - check " << nCombTrk << " comb. tracks " << endmsg;

  // loop over all combined tracks ...
  for (int i=0; i<nCombTrk;++i) {

    TrackExtended * trkExt = _allCombinedTracks[i];

    debug() << " **SelectCombinedTracks - Check Track " << trkExt << endmsg;
    
    // get the sub tracks which have been combined to form this track
    GroupTracks * group = trkExt->getGroupTracks();
    TrackExtendedVec tracks = group->getTrackExtendedVec();

    // check that there are only 2 sub tracks, as we are after Si <-> TPC mergers only
    int nTracks = int(tracks.size());

    debug() << " **SelectCombinedTracks - nTracks = " << nTracks << endmsg;
    
    if (nTracks == 2) {
      
      TrackExtended * firstTrack = tracks[0];
      TrackExtended * secondTrack = tracks[1];

      // check that the two sub tracks in question are not themselves a merger of tracks
      if ((firstTrack->getGroupTracks() == NULL) && (secondTrack->getGroupTracks() == NULL) ) {

        debug() << " **SelectCombinedTracks - firstTrack->getGroupTracks() == NULL ... "  << endmsg;

        // associate the current group to the two sub tracks  
        firstTrack->setGroupTracks(group);
        secondTrack->setGroupTracks(group);     

        // get the tracker hits ...
        TrackerHitExtendedVec firstVec = firstTrack->getTrackerHitExtendedVec();
        TrackerHitExtendedVec secondVec = secondTrack->getTrackerHitExtendedVec();

        // number of hits for first and second sub track
        int nFirst = int(firstVec.size());
        int nSecond = int(secondVec.size());

        // use these to store the min and max z positions for the combination.
        float edges[2];
        edges[0] = 1.0e+20;
        edges[1] = -1.0e+20;

        // get min and max z for the first sub track
        for (int iF=0;iF<nFirst;++iF) {
          TrackerHitExtended * trkHitExt = firstVec[iF];
	  edm4hep::ConstTrackerHit trkHit = trkHitExt->getTrackerHit();
          float zpos = float(trkHit.getPosition()[2]);
          if (zpos>edges[1])
            edges[1] = zpos;
          if (zpos<edges[0])
            edges[0] = zpos;      
        }

        // get min and max z for the second sub track
        for (int iS=0;iS<nSecond;++iS) {
          TrackerHitExtended * trkHitExt = secondVec[iS];
	  edm4hep::ConstTrackerHit trkHit = trkHitExt->getTrackerHit();
          float zpos = float(trkHit.getPosition()[2]);
          if (zpos>edges[1])
            edges[1] = zpos;
          if (zpos<edges[0])
            edges[0] = zpos;
        }

        // record the z extreams to the group ...
        group->setEdges(edges);
        // ... and add the combined track to the list.
        _trkImplVec.push_back(trkExt);
        debug() << " add track " << trkExt << " to combined final list " << endmsg;
        
	if (_debug >= 2) {
          int iopt = 1;
	  // here it is assumed that the tpc tracks is the secondTrack ...
          // PrintOutMerging(secondTrack,firstTrack,iopt);
        }       
      }
    }
    else { // if(nTracks>2) 
      debug() << " *****************  SelectCombinedTracks: MORE THAN TWO TRACKS " << nCombTrk << endmsg;
    }
  }
}

void FullLDCTrackingAlg::AddNotCombinedTracks() {  
  
  int nTPCTrk = int(_allTPCTracks.size());
  int nSiTrk = int(_allSiTracks.size());
  
  // we need some buffer vector
  TrackExtendedVec allMergedTracks;
  allMergedTracks.clear();
  
  // forcing merging of Si and TPC track segments
  if (_forceMerging==1) { 

    // loop over all TPC tracks
    for (int i=0;i<nTPCTrk;++i) {

      // get the tracks assigned to this TPC track
      TrackExtended * trkExtTPC = _allTPCTracks[i];
      GroupTracks * groupTPC = trkExtTPC->getGroupTracks();

      // if no tracks have been grouped with this TPC track
      if (groupTPC == NULL) {

        float diffMin = 1.0e+20;

        TrackExtended * siTrkToAttach = NULL;

        // loop over all Silicon Tracks
        for (int j=0;j<nSiTrk;++j) {
          TrackExtended * trkExtSi = _allSiTracks[j];
          GroupTracks * groupSi = trkExtSi->getGroupTracks();

          // only consider ungrouped Silicon Tracks
          if (groupSi == NULL) {

            int iComp = 0;
            //      float deltaP = CompareTrk(trkExtSi,trkExtTPC,_d0CutForForcedMerging,_z0CutForForcedMerging,iComp);
            float angle(0.);
            float angleSignificance(0.);
            
            // try to merge tracks using looser cuts
            float dOmega = CompareTrkII(trkExtSi,trkExtTPC,_d0CutForForcedMerging,_z0CutForForcedMerging,iComp,angle);

            float significance = CompareTrkIII(trkExtSi,trkExtTPC,_d0CutForForcedMerging,_z0CutForForcedMerging,iComp,angleSignificance);

            //      if (deltaP < _dPCutForForcedMerging) {

            if ( ((dOmega<_dOmegaForForcedMerging) && (angle<_angleForForcedMerging)) ||
                ((significance<5)                 && (angleSignificance<5))
                ) {

              float chi2O = dOmega/_dOmegaForForcedMerging;
              float chi2A = angle/_angleForForcedMerging;
              float deltaP = chi2O*chi2O + chi2A*chi2A; 

              // if this is the best match (diffMin) set the possible merger 
              if (deltaP<diffMin) {
                diffMin = deltaP;
                siTrkToAttach = trkExtSi;
              }
            } else {
              if (_debug >= 3) {
                int  iopt = 7;
                debug() << significance << " " << angleSignificance << endmsg;
                // PrintOutMerging(trkExtTPC,trkExtSi,iopt);
              }
            }
          }
        }
        
        if (siTrkToAttach!=NULL) {

          TrackExtended * trkExtSi = siTrkToAttach;
          TrackExtended * OutputTrack = new TrackExtended();
          GroupTracks * group = new GroupTracks();

          group->addTrackExtended(trkExtSi);
          group->addTrackExtended(trkExtTPC);

          OutputTrack->setGroupTracks(group);
          //        trkExtSi->setGroupTracks(group);
          //        trkExtTPC->setGroupTracks(group);         

          OutputTrack->setOmega(trkExtTPC->getOmega());
          OutputTrack->setTanLambda(trkExtSi->getTanLambda());
          OutputTrack->setPhi(trkExtSi->getPhi());
          OutputTrack->setZ0(trkExtSi->getZ0());
          OutputTrack->setD0(trkExtSi->getD0());

          float covMatTPC[15];
          float covMatSi[15];
          float covMat[15];
          for (int iCov=0;iCov<15;++iCov) {
            covMatTPC[iCov] = trkExtTPC->getCovMatrix()[iCov];
            covMatSi[iCov] = trkExtSi->getCovMatrix()[iCov];                
            covMat[iCov] = covMatSi[iCov];
          }

          float scaling = sqrt(covMatTPC[5]/covMatSi[5]);
          covMat[5]  = covMatTPC[5];
          covMat[3]  = scaling*covMatSi[3];
          covMat[4]  = scaling*covMatSi[4];
          covMat[8]  = scaling*covMatSi[8];
          covMat[12] = scaling*covMatSi[12];

          OutputTrack->setCovMatrix(covMat);
          TrackerHitExtendedVec tpcHitVec = trkExtTPC->getTrackerHitExtendedVec();
          TrackerHitExtendedVec siHitVec = trkExtSi->getTrackerHitExtendedVec();              

          int nTPCHits = int( tpcHitVec.size());
          int nSiHits = int( siHitVec.size());        

          float edges[2];
          edges[0] = 1.0e+20;
          edges[1] = -1.0e+20;

          // find the max and min z extents from hits
          for (int iH=0;iH<nSiHits;++iH) {
            TrackerHitExtended * hitExt = siHitVec[iH];
            OutputTrack->addTrackerHitExtended(hitExt);
            hitExt->setUsedInFit(true);
	    edm4hep::ConstTrackerHit hit = hitExt->getTrackerHit();
            float zpos = float(hit.getPosition()[2]);
            if (zpos<edges[0])
              edges[0] = zpos;
            if (zpos>edges[1])
              edges[1] = zpos;
          }       
          for (int iH=0;iH<nTPCHits;++iH) {
            TrackerHitExtended * hitExt = tpcHitVec[iH];
            OutputTrack->addTrackerHitExtended(hitExt);
            hitExt->setUsedInFit(true); 
	    edm4hep::ConstTrackerHit hit = hitExt->getTrackerHit();
            float zpos = float(hit.getPosition()[2]);
            if (zpos<edges[0])
              edges[0] = zpos;
            if (zpos>edges[1])
              edges[1] = zpos;
          }

          group->setEdges(edges);
          OutputTrack->setChi2(diffMin); // will be replaced if necessary
          OutputTrack->setNDF(int(1));   // will be replaced if necessary

          _allCombinedTracks.push_back( OutputTrack );
          allMergedTracks.push_back( OutputTrack );

        }
      }
    
    } // end of loop over TPC tracks
    
    
    // check that there are some merged tracks to process
    int nMerged = int(allMergedTracks.size());

    if (nMerged>0) {
    
      // sort all merged tracks by Chi2/NDF
      // although due to the fact that NDF for these tracks is set to the value 1 above,
      // it is really only a sort on Chi2 which was really the weighted difference in angle and omega

      Sorting(allMergedTracks);

      // loop over all merged tracks
      for (int iM=0;iM<nMerged;++iM) {

        
        TrackExtended * mergedTrack = allMergedTracks[iM];
        GroupTracks * grpTrk = mergedTrack->getGroupTracks();
        
        TrackExtendedVec trkVec = grpTrk->getTrackExtendedVec();
        TrackExtended * trkTPC = NULL;
        TrackExtended * trkSi = NULL;

        int nT = int(trkVec.size());

        // only consider tracks which have been composed of excactly 2 sub tracks
        if (nT==2) {

          trkTPC = trkVec[0];
          trkSi = trkVec[1];

          GroupTracks * groupTPC = trkTPC->getGroupTracks();
          GroupTracks * groupSi  = trkSi->getGroupTracks();

          // check that both the TPC and SI track have not already been combined with other tracks ...
          if (groupTPC == NULL && groupSi == NULL) {

            // set the grouping, meaning that these tracks will not be considered further
            trkTPC->setGroupTracks( grpTrk );
            trkSi->setGroupTracks( grpTrk );

            TrackerHitExtendedVec hitVec = mergedTrack->getTrackerHitExtendedVec();
            
            int nhits = int(hitVec.size());

            int totNdf = 2*nhits - 5;
            float totChi2 = trkTPC->getChi2() + trkSi->getChi2();

            mergedTrack->setNDF( totNdf );
            mergedTrack->setChi2( totChi2 );

            if (_debug >= 2) {
              int iopt = 2;
              // PrintOutMerging(trkTPC,trkSi,iopt);
            }
            _trkImplVec.push_back( mergedTrack );
          }
        }
      }
    }
  } // end of _forceMerging
  
  
  
  // clear buffer vector
  allMergedTracks.clear();
  
  // merging split up TPC segments
  if (_mergeTPCSegments) {

    std::vector<GroupTracks*> TPCSegments;
    TPCSegments.clear();

    int nNonAssignedTPCSeg = 0;

    // loop over all TPC Tracks
    for (int i=0;i<nTPCTrk;++i) {

      TrackExtended * trkExt = _allTPCTracks[i];
      GroupTracks * group = trkExt->getGroupTracks();

      debug() << " *****************  AddNotCombinedTracks: Check track " << trkExt << " id = " << trkExt->getTrack().id()  << endmsg;
      
      // only consider those tracks which have not yet been combined
      if (group == NULL) {

        // find the min and max z extents using the hits
        TrackerHitExtendedVec currentVec = trkExt->getTrackerHitExtendedVec();
        int nCur = int(currentVec.size());
        float zmin = 1e+20;
        float zmax = -1e+20;

        for (int iCur=0;iCur<nCur;++iCur) {
          TrackerHitExtended * curTrkHitExt = currentVec[iCur];
	  edm4hep::ConstTrackerHit curTrkHit = curTrkHitExt->getTrackerHit();
          float zpos = float(curTrkHit.getPosition()[2]);
          if (zpos < zmin)
            zmin = zpos;
          if (zpos > zmax)
            zmax = zpos;
        }

        
        nNonAssignedTPCSeg++;

        // current number of TPC segment groupings
        int nGroups = int(TPCSegments.size());

        float dPtMin = 1.0e+10;
        GroupTracks * groupToAttach = NULL;
        TrackExtended * trkToAttach = NULL;

        // loop over the current TPC segment groupings
        for (int iG=0;iG<nGroups;++iG) {

          GroupTracks * segments = TPCSegments[iG];
          TrackExtendedVec segVec = segments->getTrackExtendedVec();

          // number of segments with the candidate group
          int nTrk = int(segVec.size());
          bool consider = true;

          
          if (_forbidOverlapInZTPC==1) { // if overlap in Z of the two segments is forbidden

            // loop over all tracks in the current grouping
            for (int iTrk=0;iTrk<nTrk;++iTrk) {

              TrackExtended * trkInGroup = segVec[iTrk];

              // get the number of hits from the track
              TrackerHitExtendedVec hitInGroupVec = trkInGroup->getTrackerHitExtendedVec();
              int nHitsInGrp = int(hitInGroupVec.size());

              // loop over the hits and make sure that there is no overlap in z
              for (int iHitInGrp=0;iHitInGrp<nHitsInGrp;iHitInGrp++) {

                TrackerHitExtended * xTrkExt = hitInGroupVec[iHitInGrp];
		edm4hep::ConstTrackerHit xTrk = xTrkExt->getTrackerHit();

                float xZ = float(xTrk.getPosition()[2]);
                if (xZ>zmin&&xZ<zmax) {
                  consider = false;
                  break;
                }
              }
              if (!consider)
                break; // if the candiate track's min and max z are within that of the group
            }
          }

          if (consider) {

            // again loop over the tracks in the current group
            for (int iTrk=0;iTrk<nTrk;++iTrk) {

              TrackExtended * trkInGroup = segVec[iTrk];
              int iComp = 1;
              
              // compare the tracks ... 
              float dPt = CompareTrk(trkExt,trkInGroup,_d0CutToMergeTPC,_z0CutToMergeTPC,iComp);

              // check that this tracks has the lowest delta pt  and vetomerge (i.e. fullfill momentum cut and makes sure a large number of hits have not been lost in the merger)
              if (dPt < dPtMin && !VetoMerge(trkExt,trkInGroup)) {

                dPtMin = dPt;
                groupToAttach = segments;
                trkToAttach = trkInGroup;
                if (_debug>=3) {
                  int iopt = 5;
                  // PrintOutMerging(trkExt,trkInGroup,iopt);
                }
              }
              else {
                if (_debug>=3) {
                  int iopt = 9;
                  // PrintOutMerging(trkExt,trkInGroup,iopt);
                }
              }
            }
          }
          else {
            if (_debug >= 3) {
              int iopt = 9;
              for (int iTrk=0;iTrk<nTrk;++iTrk) {
                TrackExtended * trkInGroup = segVec[iTrk];
                int iComp = 1;
                float dPt = CompareTrk(trkExt,trkInGroup,_d0CutToMergeTPC,_z0CutToMergeTPC,iComp);              
                if (dPt >= dPtMin) {          
                  // PrintOutMerging(trkExt,trkInGroup,iopt);
                }
              }
            }
          }
        }

        // check the pt cut for merging and that a group has been found to match ..
        if (dPtMin < _dPCutToMergeTPC && groupToAttach != NULL) {
          
          // add the track to the group
          groupToAttach->addTrackExtended(trkExt);
          trkExt->setGroupTracks(groupToAttach);

          // set the min and max z extents
          float zminGroup = groupToAttach->getEdges()[0];
          float zmaxGroup = groupToAttach->getEdges()[1];

          float edges[2];
          edges[0] = zmin;
          if (zminGroup<zmin)
            edges[0] = zminGroup;
          edges[1] = zmax;
          if (zmaxGroup>zmax)
            edges[1] = zmaxGroup;
          groupToAttach->setEdges(edges);
          if (_debug>=3) {
            int iopt = 3;
            // PrintOutMerging(trkExt,trkToAttach,iopt);
          }
        } else {
          
          // create a new group of segments 
          GroupTracks * newSegment = new GroupTracks(trkExt);
          trkExt->setGroupTracks(newSegment);
          debug() << " *****************  AddNotCombinedTracks: Create new TPC Segment Group for track " << trkExt << " id = " << trkExt->getTrack().id()  << endmsg;
          
          TPCSegments.push_back(newSegment);
          float edges[2];
          edges[0] = zmin;
          edges[1] = zmax;
          newSegment->setEdges(edges);
        }
      }
    }

    // At this stage all tpc segements will have been grouped.
    
    // Now try to combine the groups of TPC segments with the
    // reconstructed tracks which have already been combined into full tracks
    // containing both Si and TPC segments
    
    int nCombTrk = int(_trkImplVec.size());
    int nSegments = int(TPCSegments.size());

    //    std::cout << "Combined tracks = " << nCombTrk << endmsg;
    //    std::cout << "nSegments = " << nSegments << endmsg;

    // loop over all the TPC segment collections
    for (int iS=0;iS<nSegments;++iS) {

      GroupTracks * segments = TPCSegments[iS];
      TrackExtendedVec segVec = segments->getTrackExtendedVec();

      float zminTPCSeg = segments->getEdges()[0];
      float zmaxTPCSeg = segments->getEdges()[1];

      int nTrk = int(segVec.size());
      TrackExtended * CombTrkToAttach = NULL;
      TrackExtended * keyTrack = NULL;

      float deltaPtMin = _dPCutToMergeTPC;

      // search over the combined (good) tracks
      for (int iCTrk=0;iCTrk<nCombTrk;++iCTrk) {

        TrackExtended * combTrk = _trkImplVec[iCTrk];
        GroupTracks * groupComb = combTrk->getGroupTracks();

        bool consider = true;

        if (_forbidOverlapInZComb==1) { // if overlap in Z of the two segments is forbidden
          float zminComb = groupComb->getEdges()[0];
          float zmaxComb = groupComb->getEdges()[1];
          consider = (zminTPCSeg>zmaxComb) || (zmaxTPCSeg<zminComb);
        }
        
        // if there are not overlaps in z, if _forbidOverlapInZComb is set above
        if (consider) {

          // loop over the TPC segments in the group
          for (int iTrk=0;iTrk<nTrk;++iTrk) {

            TrackExtended * trk = segVec[iTrk];
            int iopt = 0;

            // test for compatibility
            float dPt = CompareTrk(trk,combTrk,_d0CutToMergeTPC,_z0CutToMergeTPC,iopt);
            float angleSignificance(0.);
            float significance = CompareTrkIII(trk,combTrk,_d0CutToMergeTPC,_z0CutToMergeTPC,iopt,angleSignificance);
 
            // check if this is a better match than any before
            if ( (dPt<deltaPtMin || significance <5 ) ) {
              if(VetoMerge(trk,combTrk)==false){

                // asign the track to be attached
                CombTrkToAttach = combTrk;
                keyTrack = trk;
                deltaPtMin = dPt;
              }
            }
            else {
              if (_debug>=3) {
                GroupTracks * groupCur = combTrk->getGroupTracks();
                TrackExtended * dummySi = groupCur->getTrackExtendedVec()[0];
                int iopt_temp = 8;
                // PrintOutMerging(trk,dummySi,iopt_temp);
              }
            }
          }
        }
        else {
          if (_debug>=3) {
            for (int iTrk=0;iTrk<nTrk;++iTrk) {
              TrackExtended * trk = segVec[iTrk];
              int iopt = 0;
              float dPt = CompareTrk(trk,combTrk,_d0CutToMergeTPC,_z0CutToMergeTPC,iopt);
              if (dPt>deltaPtMin) {
                GroupTracks * groupCur = combTrk->getGroupTracks();
                TrackExtended * dummySi = groupCur->getTrackExtendedVec()[0];
                int iopt_temp = 8;
                // PrintOutMerging(trk,dummySi,iopt_temp);
              }
            }
          }
        }
      }
      
      if (CombTrkToAttach != NULL) { // attach TPC segment to existing Comb Track
        GroupTracks * groupToAttach = CombTrkToAttach->getGroupTracks();          
        TrackExtended * SiCombTrk = groupToAttach->getTrackExtendedVec()[0];
        TrackExtended * TpcCombTrk = groupToAttach->getTrackExtendedVec()[1];

        if (_debug>=3) {
          int iopt = 4;
          // PrintOutMerging(keyTrack,SiCombTrk,iopt);
          iopt = 5;
          // PrintOutMerging(keyTrack,TpcCombTrk,iopt);      
        }

        for (int iTrk=0;iTrk<nTrk;iTrk++) {

          TrackExtended * segmentTrack = segVec[iTrk];
          groupToAttach->addTrackExtended( segmentTrack );
          segmentTrack->setGroupTracks( groupToAttach );

          TrackerHitExtendedVec hitVec = segmentTrack->getTrackerHitExtendedVec();

          int nHitSeg = int(hitVec.size());

          for (int iHS=0;iHS<nHitSeg;++iHS) {

            // take the hit from the segment and attach it to CombTrkToAttach
            // flagging it as not to be used in the fit
            TrackerHitExtended * hitExt = hitVec[iHS];
            hitExt->setUsedInFit(false);
            CombTrkToAttach->addTrackerHitExtended( hitExt );

          }
        }
      }
      else {
        if (nTrk==1) { // create a new group 
          GroupTracks * newGrp = new GroupTracks();
          segVec[0]->setGroupTracks(newGrp);
          newGrp->addTrackExtended(segVec[0]);
          TrackerHitExtendedVec TpcHitVec = segVec[0]->getTrackerHitExtendedVec();
          int nTpcH = int(TpcHitVec.size());
          for (int iTpcH=0;iTpcH<nTpcH;++iTpcH) {
            TpcHitVec[iTpcH]->setUsedInFit( true );
          }
          _trkImplVec.push_back(segVec[0]);
          _allNonCombinedTPCTracks.push_back(segVec[0]);
        }
        else { // several segments
          
          float zMin = 1.0e+20;
          TrackExtended * chosenTrack = NULL;

          // loop over the segments
          for (int iTrk=0;iTrk<nTrk;++iTrk) {

            TrackExtended * segment = segVec[iTrk];

            // get the lcio track which is behind this segemnt
	    edm4hep::ConstTrack track = segment->getTrack();
            ConstTrackerHitVec hitVec(track.trackerHits_begin(), track.trackerHits_end());

            debug() << "Group of orphaned TPC tracks: trying track " << track.id() << endmsg;
            
            int nHits = int(hitVec.size());

            // loop over it's hits
            for (int iH=0;iH<nHits;++iH) {
              
              float zPosi = fabs(hitVec[iH].getPosition()[2]);

              // if this segment has the hit closest to the IP so far
              if (zPosi<zMin) {
                // take this as the chosen track and break
                chosenTrack = segment;
                zMin = zPosi;
                break;
              }
            }
          }
          
          if (chosenTrack!=NULL) { // can't really ever be null.
	    debug() << "Group of orphaned TPC tracks: chosen track taken as " << chosenTrack->getTrack().id() << endmsg;
            
            // create a new group of tracks
            GroupTracks * newGroup = new GroupTracks();

            // first add the chosen track
            chosenTrack->setGroupTracks( newGroup );
            newGroup->addTrackExtended( chosenTrack );
            
            // loop over the segments ...
            for (int iTrk=0;iTrk<nTrk;++iTrk) {

              TrackExtended * segment = segVec[iTrk];
              
              TrackerHitExtendedVec hitVecS = segment->getTrackerHitExtendedVec();
              int nHitS = int(hitVecS.size());                  

              // loop over the hits for the current segment
              for (int iH=0;iH<nHitS;++iH) {
                TrackerHitExtended * trkHitExt = hitVecS[iH];

                if (segment!=chosenTrack) { // ... then don't add the hits to the fit
                  // set the relation between group and track
                  segment->setGroupTracks( newGroup );
                  newGroup->addTrackExtended( segment );
                  trkHitExt->setUsedInFit( false );
                  chosenTrack->addTrackerHitExtended( trkHitExt );                              
                }
                else {
                  trkHitExt->setUsedInFit( true );
                }
              }
            }
            _allNonCombinedTPCTracks.push_back(chosenTrack);
            _trkImplVec.push_back(chosenTrack);
          }
        }
      }
    }
    for (int iS=0;iS<nSegments;++iS) {
      GroupTracks * segments = TPCSegments[iS];
      delete segments;
    }
    TPCSegments.clear();
  }
  else { // adding all TPC segments to the list of tracks (track splitting is allowed)
    for (int i=0;i<nTPCTrk;++i) {
      TrackExtended * trkExt = _allTPCTracks[i];
      edm4hep::ConstTrack track = trkExt->getTrack();
      GroupTracks * group = trkExt->getGroupTracks();

      if (group == NULL) {
        debug() << " *****************  AddNotCombinedTracks: _mergeTPCSegments = " << _mergeTPCSegments << " : Add non combined TPC track " << trkExt << " id = " << track.id()  << endmsg;

        TrackerHitExtendedVec hitVec = trkExt->getTrackerHitExtendedVec();
        int nHTPC = int(hitVec.size());

        for (int iHTPC=0;iHTPC<nHTPC;++iHTPC) {
          hitVec[iHTPC]->setUsedInFit(true);
        }
        _trkImplVec.push_back(trkExt);
        _allNonCombinedTPCTracks.push_back( trkExt );

        GroupTracks * newGrp = new GroupTracks();
        newGrp->addTrackExtended( trkExt );
        trkExt->setGroupTracks( newGrp );

      }
    }    
  }
  
  for (int i=0;i<nSiTrk;++i) { // adding left-over Si segments to the list of tracks
    TrackExtended * trkExt = _allSiTracks[i];
    edm4hep::ConstTrack track = trkExt->getTrack();
    GroupTracks * group = trkExt->getGroupTracks();

    if (group == NULL) {
      debug() << " *****************  AddNotCombinedTracks: Add non combined Silicon Track : " << trkExt << " id = " << track.id()  << endmsg;
      
      TrackerHitExtendedVec hitVec = trkExt->getTrackerHitExtendedVec();
      int nHSi = int(hitVec.size());

      for (int iHSi=0;iHSi<nHSi;++iHSi) {
        hitVec[iHSi]->setUsedInFit(true);
      }

      _trkImplVec.push_back(trkExt);
      _allNonCombinedSiTracks.push_back( trkExt );

      GroupTracks * newGrp = new GroupTracks();
      newGrp->addTrackExtended( trkExt );
      trkExt->setGroupTracks( newGrp );   

    }
  }
  
}

void FullLDCTrackingAlg::CheckTracks() {  
  
  for(unsigned int i = 0; i< _trkImplVec.size();i++){
    TrackExtended *first = _trkImplVec[i];
    if(first==NULL)continue;
    float d0First = first->getD0();
    float z0First = first->getZ0();
    float omegaFirst = first->getOmega();
    float tanLFirst = first->getTanLambda();
    float phiFirst = first->getPhi();
    HelixClass helixFirst;
    helixFirst.Initialize_Canonical(phiFirst,d0First,z0First,omegaFirst,tanLFirst,_bField);
    float momFirst[3];
    momFirst[0]= helixFirst.getMomentum()[0];
    momFirst[1]= helixFirst.getMomentum()[1];
    momFirst[2]= helixFirst.getMomentum()[2];
    float pFirst    = sqrt(momFirst[0]*momFirst[0]+momFirst[1]*momFirst[1]+momFirst[2]*momFirst[2]);
    if(std::isnan(pFirst))continue;
    TrackerHitExtendedVec firstHitVec  = first->getTrackerHitExtendedVec();
    if(firstHitVec.size()<1)continue;
    
    for(unsigned int j = i+1; j<_trkImplVec.size();j++){
      TrackExtended *second = _trkImplVec[j];
      if(second==NULL)continue;
      float d0Second = second->getD0();
      float z0Second = second->getZ0();
      float omegaSecond = second->getOmega();
      float tanLSecond = second->getTanLambda();
      float phiSecond = second->getPhi();
      HelixClass helixSecond;
      helixSecond.Initialize_Canonical(phiSecond,d0Second,z0Second,omegaSecond,tanLSecond,_bField);
      float momSecond[3];
      momSecond[0] = helixSecond.getMomentum()[0];
      momSecond[1] = helixSecond.getMomentum()[1];
      momSecond[2] = helixSecond.getMomentum()[2];
      float pSecond    = sqrt(momSecond[0]*momSecond[0]+momSecond[1]*momSecond[1]+momSecond[2]*momSecond[2]);
      if(std::isnan(pSecond))continue;
      TrackerHitExtendedVec secondHitVec  = second->getTrackerHitExtendedVec();
      if(secondHitVec.size()<1)continue;
      if(firstHitVec.size()+secondHitVec.size()<10)continue;
      
      
      float pdot = (momFirst[0]*momSecond[0]+momFirst[1]*momSecond[1]+momFirst[2]*momSecond[2])/pFirst/pSecond;
      if(pdot<0.999)continue;
      // const float sigmaPOverPFirst  = sqrt(first->getCovMatrix()[5])/fabs(omegaFirst);
      // const float sigmaPOverPSecond = sqrt(second->getCovMatrix()[5])/fabs(omegaSecond);
      // const float deltaP = fabs(pFirst-pSecond);
      // const float sigmaDeltaP = sqrt(pFirst*sigmaPOverPFirst*pFirst*sigmaPOverPFirst+pSecond*sigmaPOverPSecond*pSecond*sigmaPOverPSecond);
      //      const float significance = deltaP/sigmaDeltaP;
      
      TrackExtended * combinedTrack = CombineTracks(first,second, _maxAllowedPercentageOfOutliersForTrackCombination, true);
      if(combinedTrack != NULL){
        const int minHits = std::min(firstHitVec.size(),secondHitVec.size());
        const int maxHits = std::max(firstHitVec.size(),secondHitVec.size());
        
        if( combinedTrack->getNDF() <= 2*maxHits+minHits-5){
          delete combinedTrack->getGroupTracks();
          delete combinedTrack;
          continue;
        }
        
        float d0    = combinedTrack->getD0();
        float z0    = combinedTrack->getZ0();
        float omega = combinedTrack->getOmega();
        float tanL  = combinedTrack->getTanLambda();
        float phi   = combinedTrack->getPhi();
        
        HelixClass helix;
        helix.Initialize_Canonical(phi,d0,z0,omega,tanL,_bField);
        float mom[3];
        mom[0]  = helix.getMomentum()[0];
        mom[1]  = helix.getMomentum()[1];
        mom[2]  = helix.getMomentum()[2];
        // float p = sqrt(mom[0]*mom[0]+mom[1]*mom[1]+mom[2]*mom[2]);
        // float chi2Sig =  (combinedTrack->getChi2() - combinedTrack->getNDF());
        // chi2Sig = chi2Sig/sqrt(combinedTrack->getNDF()*2);
        
        
        
        int nTpcFirst(0);
        int nUsedFirst(0);
        for(unsigned int ihit = 0;ihit<firstHitVec.size();ihit++){
          
          if( getDetectorID(firstHitVec[ihit]->getTrackerHit()) == lcio::ILDDetID::TPC) nTpcFirst++;
          
          if( firstHitVec[ihit]->getUsedInFit()==true ) nUsedFirst++;
        }
        
        
        int nTpcSecond(0);
        int nUsedSecond(0);
        for(unsigned int ihit = 0;ihit<secondHitVec.size();ihit++){
          if( getDetectorID(secondHitVec[ihit]->getTrackerHit()) == lcio::ILDDetID::TPC) ++nTpcSecond;
          if( secondHitVec[ihit]->getUsedInFit()==true ) ++nUsedSecond;
        }
        delete combinedTrack->getGroupTracks();
        delete combinedTrack;
      }
    }
  }
  
  
}


/*
 
 compare the following:
 
 i)   delta omega < 2 * _dOmegaForMerging
 ii)  delta d0    < d0Cut
 iii) delta z0    < z0Cut

 if above three cuts are passed return:
 
 the angle between the two tracks (return by reference)
   and
 the difference in omega
 
 */


float FullLDCTrackingAlg::CompareTrkII(TrackExtended * first, TrackExtended * second, 
                                              float d0Cut, float z0Cut,int iopt,float & Angle) {
  
  
  float result = 1.0e+20;
  Angle  = 1.0e+20; 
  float omegaFirst = first->getOmega();
  float omegaSecond = second->getOmega();
  float deltaOmega = fabs((omegaFirst-omegaSecond)/omegaSecond);
  if(deltaOmega> 2*_dOmegaForMerging)return result;
  
  
  float d0First = first->getD0();
  float z0First = first->getZ0();
  float d0Second = second->getD0();
  float z0Second = second->getZ0();
  
  bool isCloseInIP = (fabs(d0First-d0Second)<d0Cut);
  isCloseInIP = isCloseInIP && (fabs(z0Second-z0First)<z0Cut);
  
  
  if(!isCloseInIP)return result;
  
  float tanLFirst = first->getTanLambda();
  float phiFirst = first->getPhi();
  float tanLSecond = second->getTanLambda();
  float phiSecond = second->getPhi();
  float qFirst = PIOVER2 - atan(tanLFirst);
  float qSecond = PIOVER2 - atan(tanLSecond);
  
  Angle = (cos(phiFirst)*cos(phiSecond)+sin(phiFirst)*sin(phiSecond))*
  sin(qFirst)*sin(qSecond)+cos(qFirst)*cos(qSecond);
  Angle = acos(Angle);
  
  result = deltaOmega;
  
  return result;
  
}


/*
 
 compare the following:

 i)   significance d0    < d0Cut
 ii)  significance z0    < z0Cut
 iii) pdot > 0.999

 if above three cuts are passed return:
 
 the significance dAngle (return by reference)
 and
 the significance dOmega
 
 */

float FullLDCTrackingAlg::CompareTrkIII(TrackExtended * first, TrackExtended * second, 
                                               float d0Cut, float z0Cut,int iopt, float & AngleSignificance) {
  
  
  float result = 1.0e+20;
  
  float d0First = first->getD0();
  float z0First = first->getZ0();
  float omegaFirst = first->getOmega();
  float tanLFirst = first->getTanLambda();
  float phiFirst = first->getPhi();
  float qFirst = PIOVER2 - atan(tanLFirst);
  
  
  float d0Second = second->getD0();
  float z0Second = second->getZ0();
  float omegaSecond = second->getOmega();
  float tanLSecond = second->getTanLambda();
  float phiSecond = second->getPhi();
  float qSecond = PIOVER2 - atan(tanLSecond);
  
  
  //MB 2010 03
  float d0ErrFirst = sqrt(first->getCovMatrix()[0]);
  float z0ErrFirst = sqrt(first->getCovMatrix()[9]);
  //  float omegaErrFirst = sqrt(first->getCovMatrix()[5]);
  float phiErrFirst = sqrt(first->getCovMatrix()[2]);
  float qErrFirst = sqrt(cos(qFirst)*cos(qFirst)*first->getCovMatrix()[14]);
  //MB END
  //MB 2010 03
  float d0ErrSecond = sqrt(second->getCovMatrix()[0]);
  float z0ErrSecond = sqrt(second->getCovMatrix()[9]);
  //  float omegaErrSecond = sqrt(second->getCovMatrix()[5]);
  float phiErrSecond = sqrt(second->getCovMatrix()[2]);
  float qErrSecond = sqrt(cos(qSecond)*cos(qSecond)*second->getCovMatrix()[14]);
  //MB END
  
  
  //  bool isCloseInIP = (fabs(d0First-d0Second)<d0Cut);
  //isCloseInIP = isCloseInIP && (fabs(z0Second-z0First)<z0Cut);
  
  //MB 2010 03
  bool isCloseInIP = (fabs(d0First-d0Second)/sqrt(d0ErrFirst*d0ErrFirst+d0ErrSecond*d0ErrSecond)<d0Cut);
  isCloseInIP = isCloseInIP && (fabs(z0Second-z0First)/sqrt(z0ErrFirst*z0ErrFirst+z0ErrSecond*z0ErrSecond)<z0Cut);
  
  if (!isCloseInIP)return result;
  
  float Angle = (cos(phiFirst)*cos(phiSecond)+sin(phiFirst)*sin(phiSecond))*
  sin(qFirst)*sin(qSecond)+cos(qFirst)*cos(qSecond);
  
  
  HelixClass helixFirst;
  helixFirst.Initialize_Canonical(phiFirst,d0First,z0First,omegaFirst,tanLFirst,_bField);
  HelixClass helixSecond;
  helixSecond.Initialize_Canonical(phiSecond,d0Second,z0Second,omegaSecond,tanLSecond,_bField);
  
  float pFirst[3];
  float pSecond[3];
  float momFirst = 0;
  float momSecond = 0;
  
  for (int iC=0;iC<3;++iC) {
    pFirst[iC] = helixFirst.getMomentum()[iC];
    pSecond[iC] = helixSecond.getMomentum()[iC];
    momFirst += pFirst[iC]* pFirst[iC];
    momSecond += pSecond[iC]*pSecond[iC];
  }
  momFirst = sqrt(momFirst);
  momSecond = sqrt(momSecond);
  
  
  float pdot = (pFirst[0]*pSecond[0]+pFirst[1]*pSecond[1]+pFirst[2]*pSecond[2])/momFirst/momSecond;
  if(pdot<0.999)return result;
  
  const float sigmaPOverPFirst  = sqrt(first->getCovMatrix()[5])/fabs(omegaFirst);
  const float sigmaPOverPSecond = sqrt(second->getCovMatrix()[5])/fabs(omegaSecond);
  const float deltaP = fabs(momFirst-momSecond);
  const float sigmaPFirst = momFirst*sigmaPOverPFirst;
  const float sigmaPSecond = momSecond*sigmaPOverPSecond;
  const float sigmaDeltaP = sqrt(sigmaPFirst*sigmaPFirst+sigmaPSecond*sigmaPSecond);
  const float significance = deltaP/sigmaDeltaP;
  
  //MB 2010 03
  float errorAngle =sin(phiFirst)*sin(phiFirst)*phiErrFirst*phiErrFirst*cos(phiSecond)*cos(phiSecond)+
  sin(phiSecond)*sin(phiSecond)*phiErrSecond*phiErrSecond*cos(phiFirst)*cos(phiFirst)+
  sin(qFirst)*sin(qFirst)*qErrFirst*qErrFirst*cos(qSecond)*cos(qSecond)+
  sin(qSecond)*sin(qSecond)*qErrSecond*qErrSecond*cos(qFirst)*cos(qFirst)+
  cos(phiFirst)*cos(phiFirst)*phiErrFirst*phiErrFirst*(sin(phiSecond)*sin(qFirst)*sin(qSecond))*(sin(phiSecond)*sin(qFirst)*sin(qSecond))+
  cos(phiSecond)*cos(phiSecond)*phiErrSecond*phiErrSecond*(sin(phiFirst)*sin(qFirst)*sin(qSecond))*(sin(phiFirst)*sin(qFirst)*sin(qSecond))+
  cos(qFirst)*cos(qFirst)*qErrFirst*qErrFirst*(sin(phiFirst)*sin(phiSecond)*sin(qSecond))*(sin(phiFirst)*sin(phiSecond)*sin(qSecond))+
  cos(qSecond)*cos(qSecond)*qErrSecond*qErrSecond*(sin(phiFirst)*sin(phiSecond)*sin(qFirst))*(sin(phiFirst)*sin(phiSecond)*sin(qFirst));
  
  if(Angle<1.){
    errorAngle = sqrt(1./(1.-Angle*Angle)*errorAngle);
  }else{
    errorAngle = sqrt(errorAngle);
  }
  
  if(errorAngle<1.e-6)errorAngle=1.e-6;
  
  AngleSignificance = fabs(acos(Angle)/errorAngle);
    
  return significance;
  
}

/*
 
 compare the following:
 
 i)   delta d0    < d0Cut  and optionally (d0_1 + d0_2) < d0cut
 ii)  delta z0    < z0Cut
 
 if above two cuts are passed then:
 
 if ( (ptFirst<_PtCutToMergeTPC) && (ptSecond<_PtCutToMergeTPC) ) {
 
     then check the difference in momentum 
 
 else 
 
     check for cases where PatRec splits non-looping TPC tracks
     look for two tracks where total tpc hits are not more than total number
     of pad rows and that the hits on one track are close to the helix of the
     other track.
 
     Check that the angular and momentum difference meets the cuts for either hight or low  pt
     also check if the angle is very small between the tracks and that the significance of the difference in pt is less than 10
 
 
 .... 
 
 note currently this is only used in AddNotCombinedTracks
 
 */


float FullLDCTrackingAlg::CompareTrk(TrackExtended * first, TrackExtended * second, float d0Cut, float z0Cut,int iopt) {
  
  float result = 1.0e+20;
  
  float d0First = first->getD0();
  float z0First = first->getZ0();
  float omegaFirst = first->getOmega();
  float tanLFirst = first->getTanLambda();
  float phiFirst = first->getPhi();
  
  float d0Second = second->getD0();
  float z0Second = second->getZ0();
  float omegaSecond = second->getOmega();
  float tanLSecond = second->getTanLambda();
  float phiSecond = second->getPhi();
  
  bool isCloseInIP = (fabs(d0First-d0Second)<d0Cut);
  
  if (iopt>0) isCloseInIP = isCloseInIP || (fabs(d0First+d0Second)<d0Cut);
  
  isCloseInIP = isCloseInIP && (fabs(z0Second-z0First)<z0Cut);
  
  
  HelixClass helixFirst;
  helixFirst.Initialize_Canonical(phiFirst,d0First,z0First,omegaFirst,tanLFirst,_bField);
  HelixClass helixSecond;
  helixSecond.Initialize_Canonical(phiSecond,d0Second,z0Second,omegaSecond,tanLSecond,_bField);
  
  float pFirst[3];
  float pSecond[3];
  float dPminus[3];
  float dPplus[3];
  float momFirst = 0;
  float momSecond = 0;
  float momMinus = 0;
  float momPlus = 0;
  
  if ( isCloseInIP ) {
    
    for (int iC=0;iC<3;++iC) {

      pFirst[iC] = helixFirst.getMomentum()[iC];
      pSecond[iC] = helixSecond.getMomentum()[iC];

      momFirst += pFirst[iC]* pFirst[iC];
      momSecond += pSecond[iC]*pSecond[iC];

      dPminus[iC] = pFirst[iC] - pSecond[iC];
      dPplus[iC] = pFirst[iC] + pSecond[iC];

      momMinus += dPminus[iC]*dPminus[iC];
      momPlus += dPplus[iC]*dPplus[iC];
    }

    momFirst = sqrt(momFirst);
    momSecond = sqrt(momSecond);
    
    float ptFirst = sqrt(pFirst[0]*pFirst[0]+pFirst[1]*pFirst[1]);
    float ptSecond = sqrt(pSecond[0]*pSecond[0]+pSecond[1]*pSecond[1]);
    
    // if both track's pt are lower than _PtCutToMergeTPC
    if ( (ptFirst<_PtCutToMergeTPC) && (ptSecond<_PtCutToMergeTPC) ) {
      
      momMinus = sqrt(momMinus);
      momPlus = sqrt(momPlus);

      float nom = momMinus;

      // get the smaller difference for the nominator
      if (momPlus<nom && iopt>0)
        nom = momPlus;

      float den = momFirst;

      // get the smallest momentum for the denominator
      if (momSecond<momFirst)
        den = momSecond;
      
      result = nom/den;     
      
    }
    else {
      
      
      // check for cases where PatRec splits non-looping TPC tracks 
      // look for two tracks where total tpc hits are not more than total number
      // of pad rows and that the hits on one track are close to the helix of the
      // other track
      
      float dpOverP = 2.0*fabs(momFirst-momSecond)/(momFirst+momSecond);
      const float pdot = (pFirst[0]*pSecond[0]+pFirst[1]*pSecond[1]+pFirst[2]*pSecond[2])/momFirst/momSecond;
      const float sigmaPOverPFirst  = sqrt(first->getCovMatrix()[5])/fabs(omegaFirst);
      const float sigmaPOverPSecond = sqrt(second->getCovMatrix()[5])/fabs(omegaSecond);
      const float deltaP = fabs(momFirst-momSecond);
      const float sigmaPFirst = momFirst*sigmaPOverPFirst;
      const float sigmaPSecond = momSecond*sigmaPOverPSecond;
      const float sigmaDeltaP = sqrt(sigmaPFirst*sigmaPFirst+sigmaPSecond*sigmaPSecond);
      const float significance = deltaP/sigmaDeltaP;
      
      
      //compare angle between the two vectors (cos theta) and their momentum
      if( (pdot>_cosThetaCutHighPtMerge && dpOverP<_momDiffCutHighPtMerge) 
         || 
         (pdot>_cosThetaCutSoftHighPtMerge && dpOverP<_momDiffCutSoftHighPtMerge) 
         || (pdot > 0.9999 && significance <10) 
         ){
        
        
        int nTrkGrpFirst = 0;
        int nTrkGrpSecond = 0;
        ConstTrackerHitVec hitvecFirst;
        ConstTrackerHitVec hitvecSecond;
        GroupTracks * groupFirst = first->getGroupTracks();
        GroupTracks * groupSecond = second->getGroupTracks();

        // does the first track belong to a group ...
        // if it does then get all the hits from the group
        if(groupFirst!=NULL){
          
          TrackExtendedVec tracksInGroupFirst = groupFirst->getTrackExtendedVec();
          nTrkGrpFirst = int(tracksInGroupFirst.size());
          
          for (int iTrkGrp=0;iTrkGrp<nTrkGrpFirst;++iTrkGrp) {
            
            TrackExtended * trkGrp = tracksInGroupFirst[iTrkGrp];
            TrackerHitExtendedVec hitVec = trkGrp->getTrackerHitExtendedVec();
            
            for(unsigned int i =0; i<hitVec.size(); ++i){
              hitvecFirst.push_back(hitVec[i]->getTrackerHit());          
            }
          }
        }

        // does the second track belong to a group ...
        // if it does then get all the hits from the group
        if(groupSecond!=NULL){
          
          TrackExtendedVec tracksInGroupSecond = groupSecond->getTrackExtendedVec();
          nTrkGrpSecond = int(tracksInGroupSecond.size());
          
          for (int iTrkGrp=0;iTrkGrp<nTrkGrpSecond;++iTrkGrp) {
            TrackExtended * trkGrp = tracksInGroupSecond[iTrkGrp];
            TrackerHitExtendedVec hitVec = 
            trkGrp->getTrackerHitExtendedVec();
            
            for(unsigned int i=0;i<hitVec.size();++i){
              hitvecSecond.push_back(hitVec[i]->getTrackerHit());
            }
          }
        }
        
        
        // for non-looping tracks 
        int nhitsFirst  = (int)hitvecFirst.size();
        int nhitsSecond = (int)hitvecSecond.size();
        int ntpcFirst   = 0;
        int ntpcSecond  = 0;
        float hitxyz[3];
        float dist[3];
        float maxdistFirst=0;
        float maxdistSecond=0;
        int ncloseFirst = 0;
        int ncloseSecond = 0;
        float zminFirst = 99999;
        float zminSecond = 99999;
        float zmaxFirst = -99999;
        float zmaxSecond = -99999;

        
        for(int ih =0;ih<nhitsFirst;++ih){
          
          float x = (float) hitvecFirst[ih].getPosition()[0];
          float y = (float) hitvecFirst[ih].getPosition()[1];
          float z = (float) hitvecFirst[ih].getPosition()[2];
          
          if(fabs(z)<zminFirst) zminFirst=fabs(z);
          if(fabs(z)>zmaxFirst) zmaxFirst=fabs(z);
          
          float r = sqrt(x*x+y*y);
        
          // count the number of hits in the TPC for the first Track
          if(r>_tpc_inner_r) ntpcFirst++;
          
          hitxyz[0]=x;
          hitxyz[1]=y;
          hitxyz[2]=z;
          helixSecond.getDistanceToPoint(hitxyz, dist);
          
          // compare 3D distance between hit and extrapolation
          if(dist[2]>maxdistFirst) maxdistFirst=dist[2];

          // count the number of hits from the first track which are closer than _hitDistanceCutHighPtMerge to the second track
          if(dist[2]<_hitDistanceCutHighPtMerge) ncloseFirst++;
        }
        
        for(int ih =0;ih<nhitsSecond;++ih){
          
          float x = (float)hitvecSecond[ih].getPosition()[0];
          float y = (float)hitvecSecond[ih].getPosition()[1];
          float z = (float)hitvecSecond[ih].getPosition()[2];
          
          if(fabs(z)<zminSecond) zminSecond=fabs(z);
          if(fabs(z)>zmaxSecond) zmaxSecond=fabs(z);
          
          float r = sqrt(x*x+y*y);
          
          // count the number of hits in the TPC for the second Track
          if(r>_tpc_inner_r) ntpcSecond++;
          
          hitxyz[0]=x;
          hitxyz[1]=y;
          hitxyz[2]=z;
          helixFirst.getDistanceToPoint(hitxyz, dist);
          
          // compare 3D distance between hit and extrapolation
          if(dist[2]>maxdistSecond) maxdistSecond=dist[2];

          // count the number of hits from the second track which are closer than _hitDistanceCutHighPtMerge to the first track
          if(dist[2]<_hitDistanceCutHighPtMerge) ncloseSecond++;

        }
        
        // calculate the fraction of hits which are closer than _hitDistanceCutHighPtMerge
        float fcloseFirst  = (float)ncloseFirst/(float)nhitsFirst;
        float fcloseSecond = (float)ncloseSecond/(float)nhitsSecond;
        
        
        
        bool split = false;
        //std::cout << "Momenta = " << momFirst << " " << momSecond << endmsg;
        //std::cout << "MaxDist = " << maxdistSecond << " " << maxdistFirst << " " << _maxHitDistanceCutHighPtMerge << endmsg;
        //std::cout << "close   = " << fcloseSecond << " " << fcloseFirst << " " << _maxFractionOfOutliersCutHighPtMerge << endmsg;
        //std::cout << "ntpc    = " << ntpcFirst << " " << ntpcSecond << " " << _tpc_pad_height+10 << endmsg;
        //std::cout << "pdot    = " << pdot << " significance " << significance << endmsg;
        
        
        // SJA:FIXME: try to fit the two tracks, without checking the number of hits which are close !!!!!!!!!!!
        TrackExtended * combinedTrack = CombineTracks(first,second, _maxAllowedPercentageOfOutliersForTrackCombination, true);

        // check that no more than 5 hits have been discared in the fit of the combined track
        if(combinedTrack != NULL){
          //std::cout << "CombinedTrack " << combinedTrack->getNDF() << " c.f. " << first->getNDF()+second->getNDF()+5 << endmsg;
          if(combinedTrack->getNDF()+10>first->getNDF()+second->getNDF()+5){
            split = true;
            dpOverP = 0;
            //std::cout << " Forcing MERGE " << endmsg;
          }
          delete combinedTrack->getGroupTracks();
          delete combinedTrack;
        } 
        else {
          //std::cout << "Could not combine track " << endmsg;
          // if close in pt and the fraction of hits matching the helix extraolations is greater than _maxFractionOfOutliersCutHighPtMerge
          if(significance<5 && fcloseFirst>_maxFractionOfOutliersCutHighPtMerge){
            split = true;
            dpOverP = 0;
            //      int overlap = SegmentRadialOverlap(first,second);
            //std::cout << " Forcing MERGE " << overlap << endmsg;
          }
        }
        
        // criteria for split track
        // old criterion
        // check the maximum deviation of the hits from the helix extrapolations, and fraction of hits matching the extrapolation based on a distance cut
        if( maxdistSecond < _maxHitDistanceCutHighPtMerge && maxdistFirst < _maxHitDistanceCutHighPtMerge 
           && 
           (fcloseSecond > _maxFractionOfOutliersCutHighPtMerge || fcloseFirst > _maxFractionOfOutliersCutHighPtMerge) 
           ){
          
          split = true;        
          
        }
        
        if(split){
          result = dpOverP;
        }
        
      }
    }
  }
  
  return result;
  
}

void FullLDCTrackingAlg::AddNotAssignedHits() {
  
  
  // currently any non-combined Silicon or TPC tracks are added to the list of tracks meaning their hits will not be available to be picked up.
  // it might be preferable to drop these tracks, at least for track in Silicon with very bad Chi2 probabilities, and allow their hits to be re-used here ...

  
  
  // Creating helix extrapolations
  //  if (_assignSETHits>0||_assignETDHits>0)
  //    CreateExtrapolations();
  
  //  if (_assignSETHits>0) { // Assignment of SET Hits
  //    
  //    const gear::GearParameters& pSETDet = Global::GEAR->getGearParameters("SET");  
  //    int nLayersSET = int(pSETDet.getDoubleVals("SETLayerRadius").size());
  //    
  //    int nSETHits = _allSETHits.size();
  //    std::vector<TrackerHitExtendedVec> SETHits;
  //    SETHits.resize(nLayersSET);
  //    
  //    for (int iSET=0;iSET<nSETHits;++iSET) {
  //      TrackerHitExtended * trkHit = _allSETHits[iSET];
  //      edm4hep::ConstTrackerHit hit = trkHit->getTrackerHit();
  //      int layer = getLayerID(trkHit);
  //      if (layer>=0&&layer<nLayersSET) 
  //        SETHits[layer].push_back(trkHit);   
  //    }
  //    for (int iL=0; iL<nLayersSET; ++iL) { // loop over SET layers
  //      TrackerHitExtendedVec hitVec = SETHits[iL];
  //      int refit = 1;
  //      AssignOuterHitsToTracks(hitVec,_distCutForSETHits,refit);
  //    }
  //  }
  
  //  if (_assignETDHits>0) { // Assignment of ETD Hits
  //    
  //    const gear::GearParameters& pETDDet = Global::GEAR->getGearParameters("ETD");  
  //    int nLayersETD = int(pETDDet.getDoubleVals("ETDLayerZ").size());
  //    
  //    int nETDHits = _allETDHits.size();
  //    std::vector<TrackerHitExtendedVec> ETDHits;
  //    ETDHits.resize(nLayersETD);
  //    
  //    for (int iETD=0;iETD<nETDHits;++iETD) {
  //      TrackerHitExtended * trkHit = _allETDHits[iETD];
  //      edm4hep::ConstTrackerHit hit = trkHit->getTrackerHit();
  //      int layer = getLayerID(trkHit);
  //      if (layer>=0 && layer < nLayersETD) 
  //        ETDHits[layer].push_back(trkHit);
  //    }
  //    for (int iL=0; iL<nLayersETD; ++iL) { // loop over ETD layers
  //      TrackerHitExtendedVec hitVec = ETDHits[iL];
  //      int refit = 0;
  //      AssignOuterHitsToTracks( hitVec, _distCutForETDHits, refit );
  //    }
  //    
  //  }
  
  //  // Cleaning up helix extrapolations
  //  if (_assignSETHits>0||_assignETDHits>0)
  //    CleanUpExtrapolations();
  
  
  
  
  
   if (_assignSETHits>0) { // Assignment of SET Hits
     debug() << "Assign SET hits *********************************" << endmsg;
     
     // Creating helix extrapolations
     CreateExtrapolations();
     
     int nSETHits = _allSETHits.size();
     std::vector<TrackerHitExtendedVec> SETHits;

     SETHits.resize(_nLayersSET);
 
     for (int iSET=0;iSET<nSETHits;++iSET) {
       TrackerHitExtended * trkHitExt = _allSETHits[iSET];
       edm4hep::ConstTrackerHit trkHit = trkHitExt->getTrackerHit();
       int layer = getLayerID(trkHit);
       if (layer>=0 && (unsigned)layer < _nLayersSET)
         SETHits[layer].push_back(trkHitExt);
     }
     for (unsigned iL=0; iL< _nLayersSET; ++iL) { // loop over SET layers
       TrackerHitExtendedVec hitVec = SETHits[iL];
       int refit = 1;
       if(hitVec.empty() == false) AssignOuterHitsToTracks(hitVec,_distCutForSETHits,refit);
     }
   }

  
  CleanUpExtrapolations();
  
  
  if (_assignSITHits>0) { // Treatment of left-over SIT hits 
    debug() << "Assign SIT hits *********************************" << endmsg;
    
    std::vector<TrackerHitExtendedVec> nonAssignedSITHits;    
    nonAssignedSITHits.resize(_nLayersSIT);
    
    int nSITHits = int(_allSITHits.size());    
    
    // loop over all SIT hits ...
    for (int iH=0;iH<nSITHits;++iH) {
      
      TrackerHitExtended * trkHitExt = _allSITHits[iH];
      TrackExtended * trkExt = trkHitExt->getTrackExtended();
      
      // check if this hit has not already been assigned to a track
      if (trkExt == NULL) {
	edm4hep::ConstTrackerHit trkHit = trkHitExt->getTrackerHit();
        
        int layer = getLayerID(trkHit);
        
        if (layer >=0 && (unsigned)layer < _nLayersSIT) {
          nonAssignedSITHits[layer].push_back(trkHitExt);
        }
      }
    }       
    
    for (int iL=_nLayersSIT-1;iL>=0;--iL) { // reverse loop over layers in Si
      
      TrackerHitExtendedVec hitVec = nonAssignedSITHits[iL];
      
      if ( hitVec.empty() == false ) {
        debug() << "AddNotAssignedHits : Try to assign hits from layer " << iL << " : Number of hits = " <<  hitVec.size() << endmsg;
        AssignSiHitsToTracks(hitVec, _distCutForSITHits);
        
      }
      
    }
  }
  
  if (_assignFTDHits>0) { // Treatment of left-over FTD hits
    debug() << "Assign FTD hits *********************************" << endmsg;
    
    std::vector<TrackerHitExtendedVec> nonAssignedFTDHits;
    nonAssignedFTDHits.resize(_nLayersFTD);
    
    int nFTDHits = int(_allFTDHits.size());
    
    // loop over all FTD hits ...
    for (int iH=0;iH<nFTDHits;++iH) {
      
      TrackerHitExtended * trkHitExt = _allFTDHits[iH];
      TrackExtended * trkExt = trkHitExt->getTrackExtended();
      
      // check if this hit has not already been assigned to a track
      if (trkExt == NULL) {
	edm4hep::ConstTrackerHit trkHit = trkHitExt->getTrackerHit();
                
        // get the layer number
        int layer = getLayerID(trkHit);
        int petalIndex = getModuleID(trkHit);
        
        if ( _petalBasedFTDWithOverlaps == true ) {
          
          // as we are dealing with staggered petals we will use 2*nlayers in each directions +/- z
          // the layers will follow the even odd numbering of the petals 
          if ( petalIndex % 2 == 0 ) {
            layer = 2*layer;
          }
          else {
            layer = 2*layer + 1;
          }
          
        }
        
        if (layer >=0 && layer < (int)_nLayersFTD)
          nonAssignedFTDHits[layer].push_back(trkHitExt);
      }
    }
    for (int iL=_nLayersFTD-1;iL>=0;--iL) {
      if ( nonAssignedFTDHits[iL].empty() == false ) {
        
        TrackerHitExtendedVec hitVec = nonAssignedFTDHits[iL];
        debug() << "AddNotAssignedHits : Try to assign hits from layer " << iL << " : Number of hits = " <<  hitVec.size() << endmsg;
        AssignSiHitsToTracks(hitVec, _distCutForFTDHits);     
        
      }
    }
  }
  
  
  
  if (_assignVTXHits>0) { // Treatment of left-over VTX hits
    debug() << "Assign VXD hits *********************************" << endmsg;
    
    std::vector<TrackerHitExtendedVec> nonAssignedVTXHits;
    nonAssignedVTXHits.resize(_nLayersVTX);
    
    int nVTXHits = int(_allVTXHits.size());
    
    // loop over all VXD hits ...
    for (int iH=0;iH<nVTXHits;++iH) {
      
      TrackerHitExtended * trkHitExt = _allVTXHits[iH];
      TrackExtended * trkExt = trkHitExt->getTrackExtended();
      
      // check if this hit has not already been assigned to a track
      if (trkExt == NULL) {
	edm4hep::ConstTrackerHit trkHit = trkHitExt->getTrackerHit();
        
        int layer = getLayerID(trkHit);
        
        if (layer >=0 && layer < (int)_nLayersVTX)
          nonAssignedVTXHits[layer].push_back(trkHitExt);
      }
    }
    for (int iL=_nLayersVTX-1;iL>=0;--iL) {
      if ( nonAssignedVTXHits[iL].empty() == false ) {
	TrackerHitExtendedVec hitVec = nonAssignedVTXHits[iL];
        debug() << "AddNotAssignedHits : Try to assign hits from layer " << iL << " : Number of hits = " <<  hitVec.size() << endmsg;
        AssignSiHitsToTracks(hitVec, _distCutForVTXHits);     
      }
    }
  }
  
  debug() << "Assign TPC hits *********************************" << endmsg;
  
  if (_assignTPCHits) {// Treatment of left-over TPC hits
    TrackerHitExtendedVec nonAssignedTPCHits;
    int nTPCHits = int(_allTPCHits.size());
    for (int iH=0;iH<nTPCHits;++iH) {
      TrackerHitExtended * trkHitExt = _allTPCHits[iH];
      TrackExtended * trkExt = trkHitExt->getTrackExtended();
      if (trkExt == NULL) {
        nonAssignedTPCHits.push_back(trkHitExt);
      }
    }
    debug() << "AddNotAssignedHits : Number of Non Assigned TPC hits = " <<  nonAssignedTPCHits.size() << endmsg;
    AssignTPCHitsToTracks(nonAssignedTPCHits, _distCutForTPCHits);
  }
  
  
}


void FullLDCTrackingAlg::CreateExtrapolations() {
  
  _trackExtrapolatedHelix.clear();
  
  int nTrk = int(_trkImplVec.size());
  for (int iTrk=0; iTrk<nTrk; ++iTrk) {
    TrackExtended * trk = _trkImplVec[iTrk];
    HelixClass * helix = GetExtrapolationHelix( trk );
    _trackExtrapolatedHelix[trk] = helix;
  }
  
}

void FullLDCTrackingAlg::CleanUpExtrapolations() {
  
  int nTrk = int(_trkImplVec.size());
  for (int iTrk=0; iTrk<nTrk; ++iTrk) {
    TrackExtended * trk = _trkImplVec[iTrk];
    HelixClass * helix =  _trackExtrapolatedHelix[trk];
    delete helix;
  }  
  _trackExtrapolatedHelix.clear();
  
}

void FullLDCTrackingAlg::AssignOuterHitsToTracks(TrackerHitExtendedVec hitVec, float dcut, int refit) {

  debug() << "FullLDCTrackingAlg::AssignOuterHitsToTracks dcut = " << dcut << endmsg;
  
  // get the number of hits to try, and the number of final tracks to which the tracks will be attached
  int nHits = int(hitVec.size());
  int nTrk = int(_trkImplVec.size());
  
  // create maps to record which tracks and tracker hits are flagged for assignment
  std::map <TrackExtended*,bool> flagTrack;
  std::map <TrackerHitExtended*,bool> flagHit;

  // vector to hold the matchups and the distance of closest approach.
  TrackHitPairVec pairs;

  flagTrack.clear();
  flagHit.clear();
  pairs.clear();
  
  // loop over all hits under consideration ...
  for (int iH=0;iH<nHits;++iH) {

    float pos[3];
    
    TrackerHitExtended * trkHitExt = hitVec[iH];
    edm4hep::ConstTrackerHit hit = trkHitExt->getTrackerHit();
    
    for (int ip=0;ip<3;++ip) {
      pos[ip] = float(hit.getPosition()[ip]);
    }
            
    // loop over all tracks in same z-half ...
    for (int iT=0;iT<nTrk;++iT) {
      
      TrackExtended * trkExt = _trkImplVec[iT];
      float tanLambda = trkExt->getTanLambda();           
      float product = pos[2]*tanLambda;
      // check that the hit and track are in the same z-half, which won't work for the rare cases of something going backwards ...
      
      if (product>0) {
      
        // use the previously created trackextrapolations for the
        HelixClass * helix = _trackExtrapolatedHelix[trkExt];
        
        // skip if the extrapolations failed
        if (helix==0) {
          debug() << "helix extrapolation failed for trkExt" << endmsg;
          continue;
        }
        
        float distance = helix->getDistanceToPoint(pos,dcut);
        
        debug() << "for helix extrapolation " << helix << " distance = " << distance << endmsg;
        
        // check the distance is less than the steerable cut ...
        if (distance<dcut) {
          debug() << "for helix extrapolation " << helix << " distance = " << distance << endmsg;
          
          // ... if so create the association and flag the hit and track
          TrackHitPair * trkHitPair =
          new TrackHitPair(trkExt,trkHitExt,distance);
          pairs.push_back(trkHitPair);
          flagTrack[trkExt] = true;
          flagHit[trkHitExt] = true;

        }
      }
    }
  }
  
  int nPairs = int(pairs.size());

  debug() << "AssignOuterHitsToTracks : Number of track hit pairs to try =  " << nPairs << endmsg;
  
  if (nPairs>0) {

    // sort the pairs on distance 
    SortingTrackHitPairs(pairs);

    for (int iP=0;iP<nPairs;++iP) {

      TrackHitPair * trkHitPair = pairs[iP];
      TrackExtended * trkExt = trkHitPair->getTrackExtended();
      TrackerHitExtended * trkHitExt = 

      trkHitPair->getTrackerHitExtended();

      // check if the track or hit is still free to be combined
      if (flagTrack[trkExt] && flagHit[trkHitExt]) {

        if (refit==0) { // just set the association
          trkExt->addTrackerHitExtended( trkHitExt );
          trkHitExt->setUsedInFit( false );
          trkHitExt->setTrackExtended( trkExt );
        }

        else {

          // get all the hits already included in the track
          TrackerHitExtendedVec hitsInTrack = trkExt->getTrackerHitExtendedVec();

          int nTotH = int(hitsInTrack.size());
          int nHitsInFit = 0;

          for (int iTH=0;iTH<nTotH;++iTH) {

            // count the number of hits used in the fit
            TrackerHitExtended * hitInTrack = hitsInTrack[iTH];
            if (hitInTrack->getUsedInFit())  {
              nHitsInFit++;
            }
          }

          int iHitInFit = 0;
          
          
          // add the previously used hits from the track to the vectors
          
          ConstTrackerHitVec trkHits;
          
          for (int iHit=0;iHit<nTotH;++iHit) {
            
            TrackerHitExtended * hitInTrack = hitsInTrack[iHit];
            if (hitInTrack->getUsedInFit()) {
	      edm4hep::ConstTrackerHit hit = hitInTrack->getTrackerHit();
              iHitInFit++;
              if(hit.isAvailable()) {
                trkHits.push_back(hit);
              }
              else{
                throw EVENT::Exception( std::string("FullLDCTrackingAlg::AssignOuterHitsToTracks: TrackerHit pointer == NULL ")  ) ;
              }
            }
          }
          
          // add the hit to be attached to the vectors
	  edm4hep::ConstTrackerHit remainHit = trkHitExt->getTrackerHit();
          iHitInFit++;
          trkHits.push_back(remainHit);
          
          
          double chi2_D;
          int ndf;
          
          
          if( trkHits.size() < 3 ) return ;
          
          // sort the hits in R
          std::vector< std::pair<float, edm4hep::ConstTrackerHit> > r2_values;
          r2_values.reserve(trkHits.size());
          
          for (ConstTrackerHitVec::iterator it=trkHits.begin(); it!=trkHits.end(); ++it) {
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
                    
          debug() << "AssignOuterHitsToTracks: Start Fitting: AddHits: number of hits to fit " << trkHits.size() << endmsg;
                    
          MarlinTrk::IMarlinTrack* marlin_trk = _trksystem->createTrack();
          
	  edm4hep::TrackState pre_fit ;
	  pre_fit.D0 = trkExt->getD0();
          pre_fit.phi = trkExt->getPhi();
          pre_fit.Z0 = trkExt->getZ0();
          pre_fit.omega = trkExt->getOmega();
          pre_fit.tanLambda = trkExt->getTanLambda();
          
          float ref[3];
          ref[0]=ref[1]=ref[2]=0.0;
          
          pre_fit.referencePoint = ref;
          
          pre_fit.location = 1/*lcio::TrackState::AtIP*/;          
          // setup initial dummy covariance matrix
          std::array<float,15> covMatrix;
          
          for (unsigned icov = 0; icov<covMatrix.size(); ++icov) {
            covMatrix[icov] = 0;
          }
          
          covMatrix[0]  = ( _initialTrackError_d0    ); //sigma_d0^2
          covMatrix[2]  = ( _initialTrackError_phi0  ); //sigma_phi0^2
          covMatrix[5]  = ( _initialTrackError_omega ); //sigma_omega^2
          covMatrix[9]  = ( _initialTrackError_z0    ); //sigma_z0^2
          covMatrix[14] = ( _initialTrackError_tanL  ); //sigma_tanl^2
          
          pre_fit.covMatrix = covMatrix;
          
          int error = MarlinTrk::createFit( trkHits, marlin_trk, &pre_fit, _bField, IMarlinTrack::backward , _maxChi2PerHit );
          
          if ( error != IMarlinTrack::success ) {
	    debug() << "FullLDCTrackingAlg::AssignOuterHitsToTracks: creation of fit fails with error " << error << endmsg;
            
            delete marlin_trk ;
            continue ;
	  }
         
          std::vector<std::pair<edm4hep::ConstTrackerHit , double> > outliers ;
          marlin_trk->getOutliers(outliers);
          
          float outlier_pct = outliers.size()/float(trkHits.size()) ;
          
          debug() << "FullLDCTrackingAlg::AssignOuterHitsToTracks: percentage of outliers " << outlier_pct << endmsg;
          
          if ( outlier_pct > _maxAllowedPercentageOfOutliersForTrackCombination) {
	    debug() << "FullLDCTrackingAlg::AssignOuterHitsToTracks: percentage of outliers " << outlier_pct << " is greater than cut maximum: " << _maxAllowedPercentageOfOutliersForTrackCombination << endmsg;
            delete marlin_trk ;
            continue ;
          }
          
          edm4hep::Vector3d point(0.,0.,0.); // nominal IP
          int return_code = 0;
          
          TrackState trkState ;
          return_code = marlin_trk->propagate(point, trkState, chi2_D, ndf ) ;
          
          delete marlin_trk ;
          
          if ( error != IMarlinTrack::success ) {
	    debug() << "FullLDCTrackingAlg::AssignOuterHitsToTracks: propagate to IP fails with error " << error << endmsg;
	    continue ;
	  }
          
          if ( ndf < 0  ) {
	    debug() << "FullLDCTrackingAlg::AssignOuterHitsToTracks: Fit failed : NDF is less that zero  " << ndf << endmsg;
	    continue ;
	  }
          
          float chi2Fit = chi2_D/float(ndf);
          
          if ( chi2Fit > _chi2FitCut ) {
	    debug() << "FullLDCTrackingAlg::AssignOuterHitsToTracks: track fail Chi2 cut of " << _chi2FitCut << " chi2 of track = " <<  chi2Fit << endmsg;
	    //why not continue? FIXME, by fucd
	    return ;
	  }
                    
          // note trackAR which is of type TrackExtended, only takes fits with ref point = 0,0,0
          trkExt->setOmega(trkState.omega);
          trkExt->setTanLambda(trkState.tanLambda);
          trkExt->setPhi(trkState.phi);
          trkExt->setD0(trkState.D0);
          trkExt->setZ0(trkState.Z0);
                    
          float cov[15];
          
          for (int i = 0 ; i<15 ; ++i) {
            cov[i] = trkState.covMatrix[i];
          }
          
          trkExt->setCovMatrix(cov);
          trkExt->setChi2(chi2_D);
          trkExt->setNDF(ndf);
          
          trkExt->addTrackerHitExtended( trkHitExt );
          trkHitExt->setTrackExtended( trkExt );
          trkHitExt->setUsedInFit( true );
          flagTrack[trkExt] = false;
          flagHit[trkHitExt] = false;
                    
          debug() << "AssignOuterHitsToTracks: Hit " << trkHitExt << " successfully assigned to track " << trkExt << endmsg;
	}
      }
    }
    
    for (int iP=0;iP<nPairs;++iP) {
      TrackHitPair * trkHitPair = pairs[iP];
      delete trkHitPair;
    }
    
    pairs.clear();
  }
}

HelixClass * FullLDCTrackingAlg::GetExtrapolationHelix( TrackExtended * track) {
  
  debug() << "FullLDCTrackingAlg::GetExtrapolationHelix called for track " << track << endmsg;
  
  HelixClass * helix = 0;
  
  // get the track state at the last hit of the track with ref point of greatest abs z
  
  GroupTracks* group = track->getGroupTracks();
  
  TrackExtendedVec trk_vec = group->getTrackExtendedVec();
  
  // location means: do not has that state
  TrackState ts_at_calo;
  ts_at_calo.location = -1;
  float z_ref_max = 0;
  
  debug() << "FullLDCTrackingAlg::GetExtrapolationHelix number of tracks associated = " << trk_vec.size() << endmsg;
  
  for (unsigned itrk=0; itrk<trk_vec.size(); ++itrk) {

    edm4hep::ConstTrack trk_lcio = trk_vec[itrk]->getTrack();
    
    if (trk_lcio.isAvailable()) {
  
      // use the tracks state at the calorimeter because that will have accounted for energy loss already
      if (hasTrackStateAt(trk_lcio, 4/*lcio::TrackState::AtCalorimeter*/)) {
        
        TrackState ts_at_last_hit = getTrackStateAt(trk_lcio, 3/*lcio::TrackState::AtLastHit*/);
        float z_ref = ts_at_last_hit.referencePoint[2];

        // make sure we use the one closest to the calo face
        if ( fabs(z_ref) >  z_ref_max) {
          z_ref_max = fabs(z_ref);
          ts_at_calo = getTrackStateAt(trk_lcio, 4/*lcio::TrackState::AtCalorimeter*/);

          debug() << "FullLDCTrackingAlg::GetExtrapolationHelix set ts_at_calo with ref_z = " << z_ref << endmsg;
	}
      }
    }
  }
  
  if (ts_at_calo.location != -1) {
    
    TrackState ts_at_calo_forIP(ts_at_calo);
        
    LCIOTrackPropagators::PropagateLCIOToNewRef(ts_at_calo_forIP,0.0,0.0,0.0);
    
    ts_at_calo_forIP.location = 1/*lcio::TrackState::AtIP*/;
    
    helix = new HelixClass();
    
    helix->Initialize_Canonical(ts_at_calo_forIP.phi,
                                ts_at_calo_forIP.D0,
                                ts_at_calo_forIP.Z0,
                                ts_at_calo_forIP.omega,
                                ts_at_calo_forIP.tanLambda,
                                _bField);
    
    debug() << "FullLDCTrackingAlg::GetExtrapolationHelix helix created at IP" << endmsg;
  }
  
  return helix;
}


void FullLDCTrackingAlg::AssignTPCHitsToTracks(TrackerHitExtendedVec hitVec,
                                                      float dcut) {
  
  int nHits = int(hitVec.size());
  int nTrk = int(_trkImplVec.size());
  
  for (int iT=0;iT<nTrk;++iT) { // loop over all tracks
    TrackExtended * foundTrack = _trkImplVec[iT];
    GroupTracks * group = foundTrack->getGroupTracks();
    TrackExtendedVec tracksInGroup = group->getTrackExtendedVec();
    int nTrkGrp = int(tracksInGroup.size());
    for (int iTrkGrp=0;iTrkGrp<nTrkGrp;++iTrkGrp) {
      TrackExtended * trkGrp = tracksInGroup[iTrkGrp];
      TrackerHitExtendedVec hitVecGrp = trkGrp->getTrackerHitExtendedVec();
      int nHits_Grp = int(hitVecGrp.size());
      float zMin = 1.0e+20;
      float zMax = -1.0e+20;
      float startPoint[3] = {0.,0.,0.};
      float endPoint[3]   = {0.,0.,0.};
      for (int iH=0;iH<nHits_Grp;++iH) {
        TrackerHitExtended * trkHitExt = hitVecGrp[iH];
        float pos[3] = {0.,0.,0.};
        for (int iC=0;iC<3;++iC) 
          pos[iC] = float(trkHitExt->getTrackerHit().getPosition()[iC]);         
        if (pos[2]>zMax) {
          zMax = pos[2];
          for (int iC=0;iC<3;++iC)
            endPoint[iC] = pos[iC];           
        }
        if (pos[2]<zMin) {
          zMin = pos[2];
          for (int iC=0;iC<3;++iC)
            startPoint[iC] = pos[iC];
        }
      }
      trkGrp->setStart(startPoint);
      trkGrp->setEnd(endPoint);
    }
  }
  
  
  // replace previous version with faster loop ordering
  
  std::vector<float> minDistances(nHits, dcut);
  std::vector<TrackExtended*> tracksToAttach(nHits);
  std::vector< std::vector<float> > HitPositions(nHits);
  std::vector<int> HitSign(nHits);//Positive or Negative side
  for (int iH=0;iH<nHits;++iH) { // loop over leftover TPC hits
    tracksToAttach[iH]=NULL;
    //Get all TrackerHit positions, so we only have to get them once
    edm4hep::ConstTrackerHit temphit = hitVec[iH]->getTrackerHit();
    const edm4hep::Vector3d temppos = temphit.getPosition();
    HitPositions[iH].push_back(float(temppos[0]));
    HitPositions[iH].push_back(float(temppos[1]));
    HitPositions[iH].push_back(float(temppos[2]));
    HitSign[iH]=std::signbit(temppos[2]);
  }    
  
  debug() << "AssignTPCHitsToTracks: Starting loop " << nTrk << " tracks   and  " << nHits << " hits" << endmsg;
  
  for (int iT=0;iT<nTrk;++iT) { // loop over all tracks
    TrackExtended * foundTrack = _trkImplVec[iT];
    int tanlambdaSign = std::signbit(foundTrack->getTanLambda());//we only care about positive or negative
    GroupTracks * group = foundTrack->getGroupTracks();
    TrackExtendedVec tracksInGroup = group->getTrackExtendedVec();
    int nTrkGrp = int(tracksInGroup.size());
    
    for (int iTrkGrp=0;iTrkGrp<nTrkGrp;++iTrkGrp) {
      TrackExtended * trkGrp = tracksInGroup[iTrkGrp];
      float tanLambda = trkGrp->getTanLambda();
      float omega = trkGrp->getOmega();
      float d0 = trkGrp->getD0();
      float z0 = trkGrp->getZ0();
      float phi0 = trkGrp->getPhi();
      float startPointZ = trkGrp->getStart()[2];
      float endPointZ   = trkGrp->getEnd()[2];
      HelixClass helix;
      helix.Initialize_Canonical(phi0,d0,z0,omega,tanLambda,_bField);
      float OnePFivehalfPeriodZ = 1.5*fabs(acos(-1.)*tanLambda/omega);
      
      for (int iH=0;iH<nHits;++iH) { // loop over leftover TPC hits
        
        //check if the hit and the track or on the same side
        //xor return 1, if hits are different
        if ( tanlambdaSign^HitSign[iH] ) continue;
        
        float DeltaStart = fabs(HitPositions[iH][2]-startPointZ);
        float DeltaEnd = fabs(HitPositions[iH][2]-endPointZ);
        bool consider = DeltaStart <= OnePFivehalfPeriodZ;
        consider = consider || (DeltaEnd <= OnePFivehalfPeriodZ);
        consider = consider || ( (HitPositions[iH][2]>=startPointZ) && (HitPositions[iH][2]<=endPointZ) );
        
        if(consider){
          float distance = helix.getDistanceToPoint(HitPositions[iH], minDistances[iH]);
          if (distance < minDistances[iH]) {
            minDistances[iH] = distance;
            tracksToAttach[iH] = foundTrack;
          }
        }
      } // loop over leftover TPC hits
    } //groups in tracks
  } // loop over all tracks
  
  for (int iH=0;iH<nHits;++iH) {
    TrackerHitExtended * trkHitExt = hitVec[iH];
    if (tracksToAttach[iH]!=NULL) {
      tracksToAttach[iH]->addTrackerHitExtended(trkHitExt);
      trkHitExt->setTrackExtended( tracksToAttach[iH] );
      trkHitExt->setUsedInFit( false );
    }
  }
  
  debug() << " Fast loop done " << endmsg;
  
  
  //     for (int iH=0;iH<nHits;iH++) { // loop over leftover TPC hits
  //    TrackerHitExtended * hitExt = hitVec[iH];
  //    float pos[3];
  //    for (int ip=0;ip<3;++ip) 
  //        pos[ip] = float(hitExt->getTrackerHit().getPosition()[ip]);
  //    float minDist = 1.0e+20;
  //    TrackExtended * trackToAttach = NULL;
  //    for (int iT=0;iT<nTrk;++iT) { // loop over all tracks
  //        TrackExtended * foundTrack = _trkImplVec[iT];
  //        float tanLambdaFound = foundTrack->getTanLambda();
  //        float product = tanLambdaFound*pos[2];
  //        if (product>0) {
  //          GroupTracks * group = foundTrack->getGroupTracks();
  //          TrackExtendedVec tracksInGroup = group->getTrackExtendedVec();
  //          int nTrkGrp = int(tracksInGroup.size());
  //          for (int iTrkGrp=0;iTrkGrp<nTrkGrp;++iTrkGrp) {
  //            TrackExtended * trkGrp = tracksInGroup[iTrkGrp];
  //            float tanLambda = trkGrp->getTanLambda();
  //            float omega = trkGrp->getOmega();
  //            float d0 = trkGrp->getD0();
  //            float z0 = trkGrp->getZ0();
  //            float phi0 = trkGrp->getPhi();
  //            float dist[3];
  //            float startPoint[3];
  //            float endPoint[3];
  //            for (int iC=0;iC<3;++iC) {
  //              startPoint[iC] = trkGrp->getStart()[iC];
  //              endPoint[iC] = trkGrp->getEnd()[iC];
  //            }
  //            HelixClass helix;
  //            helix.Initialize_Canonical(phi0,d0,z0,omega,tanLambda,_bField);
  //            float halfPeriodZ = fabs(acos(-1.)*tanLambda/omega);
  //            helix.getDistanceToPoint(pos,dist);
  //            float DeltaStart = fabs(pos[2]-startPoint[2]);
  //            float DeltaEnd = fabs(pos[2]-endPoint[2]);
  //            bool consider = DeltaStart <= 1.5*halfPeriodZ;
  //            consider = consider || (DeltaEnd <= 1.5*halfPeriodZ);
  //            consider = consider || ( (pos[2]>=startPoint[2]) && (pos[2]<=endPoint[2]) );
  // //                 float ZMin = DeltaStart;
  // //                 if (DeltaEnd<ZMin)
  // //                   ZMin = DeltaEnd;
  //            if (dist[2]<dcut && consider && dist[2]<minDist) {
  //              minDist = dist[2];
  //              trackToAttach = foundTrack;
  //            }
  //          }
  //        }
  //    }
  
  //    if (trackToAttach!=NULL) {
  //   trackToAttach->addTrackerHitExtended(hitExt);
  //  hitExt->setTrackExtended( trackToAttach );
  //  hitExt->setUsedInFit( false );
  //  if(trackToAttach!=tracksToAttach[iH])std::cout << " Check Failed" << trackToAttach << "  " << tracksToAttach[iH] << endmsg;
  //
  //}
  //else {
  ///     std::cout << iH << " hit is not assigned : distance to closest track = " << minDist << endmsg;
  ///}
  //}
  //std::cout << " Slow loop done " << endmsg;
  
}

void FullLDCTrackingAlg::AssignSiHitsToTracks(TrackerHitExtendedVec hitVec,
                                                     float dcut) {
  
  int nHits = int(hitVec.size());
  int nTrk = int(_allNonCombinedTPCTracks.size());
  
  debug() << "AssignSiHitsToTracks : Number of hits to assign " <<  hitVec.size() << " : Number of available tracks = " << nTrk << endmsg;
  
  std::map <TrackExtended*,bool> flagTrack;
  std::map <TrackerHitExtended*,bool> flagHit;
  TrackHitPairVec pairs;
  flagTrack.clear();
  flagHit.clear();
  pairs.clear();
  
  for (int iH=0;iH<nHits;++iH) {
    
    float pos[3];
    TrackerHitExtended * trkHitExt = hitVec[iH];
    edm4hep::ConstTrackerHit hit = trkHitExt->getTrackerHit();
    
    for (int ip=0;ip<3;++ip) {
      pos[ip] = float(hit.getPosition()[ip]);
    }
    
    for (int iT=0;iT<nTrk;++iT) {
      
      TrackExtended * trkExt = _allNonCombinedTPCTracks[iT];
      
      float tanLambda = trkExt->getTanLambda();       
      float product = pos[2]*tanLambda;
      
      debug() << "AssignSiHitsToTracks : product =  " << product << " z hit = " << pos[2] <<  endmsg;
      
      if (product>0) {
        
        float d0 = trkExt->getD0();
        float z0 = trkExt->getZ0();
        float phi0 = trkExt->getPhi();
        float omega = trkExt->getOmega();
        tanLambda = trkExt->getTanLambda();
        
        HelixClass helix;
        helix.Initialize_Canonical(phi0,d0,z0,omega,tanLambda,_bField);
        float distance = helix.getDistanceToPoint(pos,dcut);
        
        debug() << "AssignSiHitsToTracks : distance =  " << distance << " cut = " << dcut << endmsg;
        
        if (distance<dcut) {
          TrackHitPair * trkHitPair = 
          new TrackHitPair(trkExt,trkHitExt,distance);
          pairs.push_back(trkHitPair);
          flagTrack[trkExt] = true;
          flagHit[trkHitExt] = true;
        }
      }
    }
  }
  
  int nPairs = int(pairs.size());
  
  debug() << "AssignSiHitsToTracks : Number of track hit pairs to try =  " << nPairs << endmsg;
  
  if (nPairs>0) {
    
    SortingTrackHitPairs(pairs);
    
    for (int iP=0;iP<nPairs;++iP) {
      
      TrackHitPair * trkHitPair = pairs[iP];
      TrackExtended * trkExt = trkHitPair->getTrackExtended();
      TrackerHitExtended * trkHitExt = 
      
      trkHitPair->getTrackerHitExtended();
      
      if (flagTrack[trkExt] && flagHit[trkHitExt]) {              
        
        // get the hits already assigned to the track
        TrackerHitExtendedVec hitsInTrack = trkExt->getTrackerHitExtendedVec();
        
        int nTotH = int(hitsInTrack.size());
        int nHitsInFit = 0;
        
        for (int iTH=0;iTH<nTotH;++iTH) {
          
          TrackerHitExtended * hitInTrack = hitsInTrack[iTH];
          
          // count the number of hits used in the fit
          if (hitInTrack->getUsedInFit()) {
            nHitsInFit++; 
          }
        }
        
        int iHitInFit = 0;
        
        // add the previously used hits from the track to the vectors 
        
        ConstTrackerHitVec trkHits;
        
        for (int iHit=0;iHit<nTotH;++iHit) {
          
          TrackerHitExtended * hitInTrack = hitsInTrack[iHit];
          if (hitInTrack->getUsedInFit()) {
	    edm4hep::ConstTrackerHit hit = hitInTrack->getTrackerHit();
            iHitInFit++;
            if(hit.isAvailable()) { 
              trkHits.push_back(hit);   
            }
            else{
              throw EVENT::Exception( std::string("FullLDCTrackingAlg::AssignSiHitsToTracks: TrackerHit pointer == NULL ")  ) ;
            }
          }
        }
        
        // add the hit to be attached to the vectors 
	edm4hep::ConstTrackerHit remainHit = trkHitExt->getTrackerHit();
        iHitInFit++;
        trkHits.push_back(remainHit);
        
        
        double chi2_D;
        int ndf;
        
        if( trkHits.size() < 3 ) return ;
        
        // sort the hits in R
        std::vector< std::pair<float, edm4hep::ConstTrackerHit> > r2_values;
        r2_values.reserve(trkHits.size());
        
        for (ConstTrackerHitVec::iterator it=trkHits.begin(); it!=trkHits.end(); ++it) {
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

	debug() << "AssignSiHitsToTracks: Start Fitting: AddHits: number of hits to fit " << trkHits.size() << endmsg;
        
        MarlinTrk::IMarlinTrack* marlin_trk = _trksystem->createTrack();
        
        TrackState pre_fit ;
        pre_fit.D0 = trkExt->getD0();
        pre_fit.phi = trkExt->getPhi();
        pre_fit.Z0 = trkExt->getZ0();
        pre_fit.omega = trkExt->getOmega();
        pre_fit.tanLambda = trkExt->getTanLambda();
        
        float ref[3];
        ref[0]=ref[1]=ref[2]=0.0;
        
        pre_fit.referencePoint = ref;
        
        pre_fit.location = 1/*lcio::TrackState::AtIP*/;
        
        // setup initial dummy covariance matrix
        std::array<float,15> covMatrix;

        for (unsigned icov = 0; icov<covMatrix.size(); ++icov) {
          covMatrix[icov] = 0;
        }
        
        covMatrix[0]  = ( _initialTrackError_d0    ); //sigma_d0^2
        covMatrix[2]  = ( _initialTrackError_phi0  ); //sigma_phi0^2
        covMatrix[5]  = ( _initialTrackError_omega ); //sigma_omega^2
        covMatrix[9]  = ( _initialTrackError_z0    ); //sigma_z0^2
        covMatrix[14] = ( _initialTrackError_tanL  ); //sigma_tanl^2
        
        pre_fit.covMatrix = covMatrix;
        
        int error = MarlinTrk::createFit( trkHits, marlin_trk, &pre_fit, _bField, IMarlinTrack::backward , _maxChi2PerHit );
        
        if ( error != IMarlinTrack::success ) {
          debug() << "FullLDCTrackingAlg::AssignSiHitsToTracks: creation of fit fails with error " << error << endmsg;
          
          delete marlin_trk ;
          continue ;
	}
        
        std::vector<std::pair<edm4hep::ConstTrackerHit , double> > outliers ;
        marlin_trk->getOutliers(outliers);
        
        float outlier_pct = outliers.size()/float(trkHits.size());
        
        debug()<< "FullLDCTrackingAlg::AssignSiHitsToTracks: percentage of outliers " << outlier_pct << endmsg;
        
        if ( outlier_pct > _maxAllowedPercentageOfOutliersForTrackCombination) {
	  debug() << "FullLDCTrackingAlg::AssignSiHitsToTracks: percentage of outliers " << outlier_pct << " is greater than cut maximum: " << _maxAllowedPercentageOfOutliersForTrackCombination << endmsg;
          delete marlin_trk ;
          continue ;
          
        }
                
        edm4hep::Vector3d point(0.,0.,0.); // nominal IP
        int return_code = 0;
        
        TrackState trkState ;
        return_code = marlin_trk->propagate(point, trkState, chi2_D, ndf ) ;
        
        delete marlin_trk ;
        
        if ( error != IMarlinTrack::success ) {
	  debug() << "FullLDCTrackingAlg::AssignSiHitsToTracks: propagate to IP fails with error " << error << endmsg;
	  continue ;
	}
        
        if ( ndf < 0  ) {
	  debug() << "FullLDCTrackingAlg::AssignSiHitsToTracks: Fit failed : NDF is less that zero  " << ndf << endmsg;
	  continue ;
	}
        
        float chi2Fit = chi2_D/float(ndf);
        
        if ( chi2Fit > _chi2FitCut ) {
	  debug() << "FullLDCTrackingAlg::AssignSiHitsToTracks: track fail Chi2 cut of " << _chi2FitCut << " chi2 of track = " <<  chi2Fit << endmsg;
	  continue ;
	}
                
        // note trackAR which is of type TrackExtended, only takes fits with ref point = 0,0,0 
        trkExt->setOmega(trkState.omega);
        trkExt->setTanLambda(trkState.tanLambda);
        trkExt->setPhi(trkState.phi);
        trkExt->setD0(trkState.D0);
        trkExt->setZ0(trkState.Z0);
        
        float cov[15];
        
        for (int i = 0 ; i<15 ; ++i) {
          cov[i] = trkState.covMatrix[i];
        }
        
        trkExt->setCovMatrix(cov);
        trkExt->setChi2(chi2_D);
        trkExt->setNDF(ndf);
        
        trkExt->addTrackerHitExtended( trkHitExt );
        trkHitExt->setTrackExtended( trkExt );
        trkHitExt->setUsedInFit( true );
        flagTrack[trkExt] = false;
        flagHit[trkHitExt] = false;
                
        debug() << "AssignSiHitsToTracks: Hit " << trkHitExt << " successfully assigned to track " << trkExt << endmsg;
      }
    }
    
    for (int iP=0;iP<nPairs;++iP) {
      TrackHitPair * trkHitPair = pairs[iP];
      delete trkHitPair;
    }
    
    pairs.clear();
    
  }
}

/*
void FullLDCTrackingAlg::PrintOutMerging(TrackExtended * firstTrackExt, TrackExtended * secondTrackExt, int iopt) {
  // iopt = 1 merged Si and TPC merging
  // iopt = 2 merged Si and TPC using forced merging
  // iopt = 3 merged TPC segments merging
  // iopt = 4 merged Comb Si and TPC merging
  // iopt = 5 merged Comb TPC and TPC merging
  // iopt = 6 unmerged TPC and Si segments ( when using soft merging)
  // iopt = 7 unmerged TPC and Si segments ( when using forced merging)
  // iopt = 8 unmerged Comb and TPC
  // iopt = 9 unmerged TPC segments
  
  char strg[200];
  
  try {
    
    Track * firstTrack = firstTrackExt->getTrack();
    Track * secondTrack = secondTrackExt->getTrack();
    
    std::string firstColName  = _TPCTrackMCPCollName;
    std::string secondColName = _TPCTrackMCPCollName;
    
    if (iopt==1) {
      secondColName = _SiTrackMCPCollName;
    }
    else if (iopt==2) {
      secondColName = _SiTrackMCPCollName;
    }
    else if (iopt==3) {
      secondColName = _TPCTrackMCPCollName;
    }
    else if (iopt==4) {
      secondColName = _SiTrackMCPCollName;
    }
    else if (iopt==5) {
      secondColName = _TPCTrackMCPCollName;
    }
    else if (iopt==6) {
      secondColName = _SiTrackMCPCollName;
    }
    else if (iopt==7) {
      secondColName = _SiTrackMCPCollName;
    }
    else if (iopt==8) {
      secondColName = _SiTrackMCPCollName;
    }    
    else {
      secondColName = _TPCTrackMCPCollName;
    }
    
    
    // get the relation collections, this should be done only once for each event ...
    LCCollection * firstCol = _evt->getCollection(firstColName.c_str());
    LCCollection * secondCol = _evt->getCollection(secondColName.c_str());
    
    // get navigators
    LCRelationNavigator firstNav(firstCol);
    LCRelationNavigator secondNav(secondCol);

    
    // get the MCParticles with the greatest weight for the first track
    LCObjectVec firstVec = firstNav.getRelatedToObjects(firstTrack);
    FloatVec firstWeights = firstNav.getRelatedToWeights(firstTrack);

    LCObject * firstMCP = NULL;

    float firstWght = 0;
    int nObj = firstVec.size();
    for (int iObj=0;iObj<nObj;++iObj) {
      if (firstWeights[iObj]>firstWght) {
        firstWght = firstWeights[iObj];
        firstMCP = firstVec[iObj];
      }
    }
    
    
    // get the MCParticles with the greatest weight for the second track
    LCObjectVec secondVec = secondNav.getRelatedToObjects(secondTrack);
    FloatVec secondWeights = secondNav.getRelatedToWeights(secondTrack);

    LCObject * secondMCP = NULL;

    float secondWght = 0;
    nObj = secondVec.size();
    for (int iObj=0;iObj<nObj;++iObj) {
      if (secondWeights[iObj]>secondWght) {
        secondWght = secondWeights[iObj];
        secondMCP = secondVec[iObj];
      }
    }
    

    // get the track parameters for both tracks and get the 3-momentum using the HelixClass
    float d0First = firstTrackExt->getD0();
    float z0First = firstTrackExt->getZ0();
    float omegaFirst = firstTrackExt->getOmega();
    float tanLFirst = firstTrackExt->getTanLambda();
    float phi0First = firstTrackExt->getPhi();
    
    float d0Second = secondTrackExt->getD0();
    float z0Second = secondTrackExt->getZ0();
    float omegaSecond = secondTrackExt->getOmega();
    float tanLSecond = secondTrackExt->getTanLambda();
    float phi0Second = secondTrackExt->getPhi();            
    
    HelixClass firstHelix;
    firstHelix.Initialize_Canonical(phi0First,d0First,z0First,omegaFirst,tanLFirst,_bField);
    float pxFirst = firstHelix.getMomentum()[0];
    float pyFirst = firstHelix.getMomentum()[1];
    float pzFirst = firstHelix.getMomentum()[2];            
    
    HelixClass secondHelix;
    secondHelix.Initialize_Canonical(phi0Second,d0Second,z0Second,omegaSecond,tanLSecond,_bField);
    float pxSecond = secondHelix.getMomentum()[0];
    float pySecond = secondHelix.getMomentum()[1];
    float pzSecond = secondHelix.getMomentum()[2];          

    
    // get the momentum differences
    float dPx = pxFirst + pxSecond;
    float dPy = pyFirst + pySecond;
    float dPz = pzFirst + pzSecond;
    
    float dPplus  = sqrt(dPx*dPx+dPy*dPy+dPz*dPz);
    
    dPx = pxFirst - pxSecond;
    dPy = pyFirst - pySecond;
    dPz = pzFirst - pzSecond;
    
    float dPminus = sqrt(dPx*dPx+dPy*dPy+dPz*dPz);
    
    // get momentum for each track
    float pFirst  = sqrt(pxFirst * pxFirst+ pyFirst* pyFirst+ pzFirst*pzFirst);
    float pSecond = sqrt(pxSecond*pxSecond+pySecond*pySecond+pzSecond*pzSecond);
    
    float ptFirst  = sqrt(pxFirst * pxFirst+ pyFirst* pyFirst);
    float ptSecond = sqrt(pxSecond*pxSecond+pySecond*pySecond);
    
    
    const float sigmaPOverPFirst  = sqrt( firstTrackExt->getCovMatrix()[5])/fabs(omegaFirst);
    const float sigmaPOverPSecond = sqrt(secondTrackExt->getCovMatrix()[5])/fabs(omegaSecond);

    const float sigmaPFirst  =  pFirst*sigmaPOverPFirst;
    const float sigmaPSecond = pSecond*sigmaPOverPSecond;
    
    
    float den = pFirst;

    if (pSecond<pFirst)
      den = pSecond;
    
    dPplus  = dPplus/den;
    dPminus = dPminus/den; 
    
    // now check if this was a Erroneous merger ...
    if (firstMCP!=secondMCP && iopt < 6) {
      
      if (iopt==1) {
        streamlog_out(DEBUG4) << "Erroneous merging of Si and TPC segments (iopt=1) mcp first = " << firstMCP << " mcp second = " << secondMCP << " --->" << endmsg;
      }
      else if (iopt==2) {
        streamlog_out(DEBUG4) << "Erroneous merging of Si and TPC segments (iopt=2) mcp first = " << firstMCP << " mcp second = " << secondMCP << " --->" << endmsg;
      }
      else if (iopt==3) {
        streamlog_out(DEBUG4) << "Erroneous merging of TPC segments (iopt=3) mcp first = " << firstMCP << " mcp second = " << secondMCP << " ---> " << endmsg;
      }
      else if (iopt==4) {
        streamlog_out(DEBUG4) << "Erroneous merging of combSi segment with uncombTPC segment (iopt=4) mcp first = " << firstMCP << " mcp second = " << secondMCP << " ---> " << endmsg;
      }
      else {
        streamlog_out(DEBUG4) << "Erroneous merging of combTPC segment with uncombTPC segment (iopt=5) mcp first = " << firstMCP << " mcp second = " << secondMCP << " --->" << endmsg;
      }
      
      
      streamlog_out(DEBUG4) << "    p         error      pt       D0      Z0     Px      Py      Pz      wieght" << endmsg;
      
      sprintf(strg,"%7.2f +- %7.2f   %7.2f  %7.1f %7.1f %7.2f %7.2f %7.2f  ",
              pFirst, sigmaPFirst, ptFirst, d0First,z0First,pxFirst,pyFirst,pzFirst);
      streamlog_out(DEBUG4) << strg;
      
      sprintf(strg,"  %5.3f\n",firstWght);
      streamlog_out(DEBUG4) << strg;
      
      
      sprintf(strg,"%7.2f +- %7.2f   %7.2f  %7.1f %7.1f %7.2f %7.2f %7.2f  ",
              pSecond, sigmaPSecond, ptSecond, d0Second,z0Second,pxSecond,pySecond,pzSecond);
      streamlog_out(DEBUG4) << strg;
      sprintf(strg,"  %5.3f\n",secondWght);
      streamlog_out(DEBUG4) << strg;
      
      streamlog_out(DEBUG4) << "Difference in +P = " << dPplus << "  -P = " << dPminus << endmsg;
      
      TrackExtended * combinedTrack = CombineTracks(firstTrackExt,secondTrackExt, _maxAllowedPercentageOfOutliersForTrackCombination, true);
      
      if(combinedTrack != NULL){
        streamlog_out(DEBUG4) << "CombinedTrack " << combinedTrack->getNDF() << " c.f. " << firstTrackExt->getNDF()+secondTrackExt->getNDF()+5 << endmsg;
        delete combinedTrack->getGroupTracks();
        delete combinedTrack;
      }else{
        streamlog_out(DEBUG4) << "Could not combine track " << endmsg;
      }
      streamlog_out(DEBUG4) << " Overlap = " << SegmentRadialOverlap(firstTrackExt, secondTrackExt) << " veto = " << VetoMerge(firstTrackExt, secondTrackExt) << endmsg;
      
      streamlog_out(DEBUG4) << endmsg;
      
    }
    
    // ... or if it was an incorrect TPC to TPC rejection ...
    else if (firstMCP==secondMCP && ( (iopt==8) || (iopt==9) ) ) {
      
      float deltaOmega = _dOmegaForMerging;
      float deltaAngle = _angleForMerging;
      if (iopt==8) {
        streamlog_out(DEBUG4) << "Unmerged TPC and Comb segments (iopt=8) --->" << endmsg;
      }
      else {
        streamlog_out(DEBUG4) << "Unmerged TPC segments (iopt=9) --->" << endmsg;
      }
      
      float qFirst = PIOVER2 - atan(tanLFirst);
      float qSecond = PIOVER2 - atan(tanLSecond);
      float dOmega = fabs((omegaFirst-omegaSecond)/omegaSecond);
      float angle = (cos(phi0First)*cos(phi0Second)+sin(phi0First)*sin(phi0Second))*
      sin(qFirst)*sin(qSecond)+cos(qFirst)*cos(qSecond);
      angle = acos(angle);
      
      
      streamlog_out(DEBUG4) << " dOmegaCut = " << deltaOmega
      << " AngleCut = " << deltaAngle
      << " dOmega = " << dOmega
      << " angle = " << angle << endmsg;
      
      
      streamlog_out(DEBUG4) << "    p         error      pt       D0      Z0     Px      Py      Pz      wieght" << endmsg;
      
      sprintf(strg,"%7.2f +- %7.2f   %7.2f  %7.1f %7.1f %7.2f %7.2f %7.2f  ",
              pFirst, sigmaPFirst, ptFirst, d0First,z0First,pxFirst,pyFirst,pzFirst);
      streamlog_out(DEBUG4) << strg;

      sprintf(strg,"  %5.3f\n",firstWght);
      streamlog_out(DEBUG4) << strg;

      
      sprintf(strg,"%7.2f +- %7.2f   %7.2f  %7.1f %7.1f %7.2f %7.2f %7.2f  ",
              pSecond, sigmaPSecond, ptSecond, d0Second,z0Second,pxSecond,pySecond,pzSecond);
      streamlog_out(DEBUG4) << strg;

      sprintf(strg,"  %5.3f\n",secondWght);
      streamlog_out(DEBUG4) << strg;
      
      streamlog_out(DEBUG4) << "Difference in +P = " << dPplus << "  -P = " << dPminus << endmsg;
      
      TrackExtended * combinedTrack = CombineTracks(firstTrackExt,secondTrackExt, _maxAllowedPercentageOfOutliersForTrackCombination, true);
      if(combinedTrack != NULL){
        streamlog_out(DEBUG4) << "CombinedTrack " << combinedTrack->getNDF() << " c.f. " << firstTrackExt->getNDF()+secondTrackExt->getNDF()+5 << endmsg;
        delete combinedTrack->getGroupTracks();
        delete combinedTrack;
      }else{
        streamlog_out(DEBUG4) << "Could not combine track " << endmsg;
      }
      streamlog_out(DEBUG4) << " Overlap = " << SegmentRadialOverlap(firstTrackExt, secondTrackExt) << " veto = " << VetoMerge(firstTrackExt, secondTrackExt) << endmsg;
      
      streamlog_out(DEBUG4) << endmsg;
      
    }    

    // ... or if it was an incorrect TPC to Si rejection ...
    else if (firstMCP==secondMCP && ( (iopt == 6) || (iopt == 7) ) ) {
      
      float deltaOmega = _dOmegaForMerging;
      float deltaAngle = _angleForMerging;
      
      if (iopt ==6) {
        streamlog_out(DEBUG4) << "Unmerged TPC and Si segments (iopt=6) mcp first = " << firstMCP << " mcp second = " << secondMCP << " --->" << endmsg;
      }
      else {
        streamlog_out(DEBUG4) << "Unmerged TPC and Si segments (iopt=7) mcp first = " << firstMCP << " mcp second = " << secondMCP << " --->" << endmsg;
        deltaOmega = _dOmegaForForcedMerging;
        deltaAngle = _angleForForcedMerging;
      }
      
      float qFirst = PIOVER2 - atan(tanLFirst);
      float qSecond = PIOVER2 - atan(tanLSecond);
      
      float dOmega = fabs((omegaFirst-omegaSecond)/omegaSecond);
      float angle = (cos(phi0First)*cos(phi0Second)+sin(phi0First)*sin(phi0Second))*
      sin(qFirst)*sin(qSecond)+cos(qFirst)*cos(qSecond);
      angle = acos(angle);
      
      
      streamlog_out(DEBUG4) << " dOmegaCut = " << deltaOmega
      << " AngleCut = " << deltaAngle
      << " dOmega = " << dOmega
      << " angle = " << angle << endmsg;

      streamlog_out(DEBUG4) << "    p         error      pt       D0      Z0     Px      Py      Pz      wieght" << endmsg;
      
      sprintf(strg,"%7.2f +- %7.2f   %7.2f  %7.1f %7.1f %7.2f %7.2f %7.2f  ",
              pFirst, sigmaPFirst, ptFirst, d0First,z0First,pxFirst,pyFirst,pzFirst);
      streamlog_out(DEBUG4) << strg;

      sprintf(strg,"  %5.3f\n",firstWght);
      streamlog_out(DEBUG4) << strg;

      sprintf(strg,"%7.2f +- %7.2f   %7.2f  %7.1f %7.1f %7.2f %7.2f %7.2f  ",
              pSecond, sigmaPSecond, ptSecond, d0Second,z0Second,pxSecond,pySecond,pzSecond);
      streamlog_out(DEBUG4) << strg;

      sprintf(strg,"  %5.3f\n",secondWght);
      streamlog_out(DEBUG4) << strg;
      
      streamlog_out(DEBUG4) << "Difference in +P = " << dPplus << "  -P = " << dPminus << endmsg;      
      
      TrackExtended * combinedTrack = CombineTracks(firstTrackExt,secondTrackExt, _maxAllowedPercentageOfOutliersForTrackCombination, true);
      
      if(combinedTrack != NULL){
        streamlog_out(DEBUG4) << "CombinedTrack " << combinedTrack->getNDF() << " c.f. " << firstTrackExt->getNDF()+secondTrackExt->getNDF()+5 << endmsg;
        delete combinedTrack->getGroupTracks();
        delete combinedTrack;
      }else{
        streamlog_out(DEBUG4) << "Could not combine track " << endmsg;
      }
      streamlog_out(DEBUG4) << " Overlap = " << SegmentRadialOverlap(firstTrackExt, secondTrackExt) << " veto = " << VetoMerge(firstTrackExt, secondTrackExt) << endmsg;
      streamlog_out(DEBUG4) << endmsg;      
      
    }
    // ... or if it was an correct merger ...
    else if (firstMCP==secondMCP && iopt < 6 && _debug > 3) {
      //      return;
      if (iopt==1) {
        streamlog_out(DEBUG4) << "Correctly combining Si and TPC segments (iopt=1) mcp first = " << firstMCP << " mcp second = " << secondMCP << " --->" << endmsg;
      }
      else if (iopt==2) {
        streamlog_out(DEBUG4) << "Correctly merging of Si and TPC segments (iopt=2) mcp first = " << firstMCP << " mcp second = " << secondMCP << " --->" << endmsg;
      }
      else if (iopt==3) {
        streamlog_out(DEBUG4) << "Correctly merging of TPC segments (iopt=3) mcp first = " << firstMCP << " mcp second = " << secondMCP << " ---> " << endmsg;
      }
      else if (iopt==4) {
        streamlog_out(DEBUG4) << "Correctly merging of combSi segment with uncombTPC segment (iopt=4) mcp first = " << firstMCP << " mcp second = " << secondMCP << " ---> " << endmsg;
      }
      else {
        streamlog_out(DEBUG4) << "Correctly merging of combTPC segment with uncombTPC segment (iopt=5) mcp first = " << firstMCP << " mcp second = " << secondMCP << " --->" << endmsg;
      }

      streamlog_out(DEBUG4) << "    p         error      pt       D0      Z0     Px      Py      Pz      wieght" << endmsg;
      
      sprintf(strg,"%7.2f +- %7.2f   %7.2f  %7.1f %7.1f %7.2f %7.2f %7.2f  ",
              pFirst, sigmaPFirst, ptFirst, d0First,z0First,pxFirst,pyFirst,pzFirst);
      streamlog_out(DEBUG4) << strg;

      sprintf(strg,"  %5.3f\n",firstWght);
      streamlog_out(DEBUG4) << strg;
      
      sprintf(strg,"%7.2f +- %7.2f   %7.2f  %7.1f %7.1f %7.2f %7.2f %7.2f  ",
              pSecond, sigmaPSecond, ptSecond, d0Second,z0Second,pxSecond,pySecond,pzSecond);
      streamlog_out(DEBUG4) << strg;

      sprintf(strg,"  %5.3f\n",secondWght);
      streamlog_out(DEBUG4) << strg;
      
      streamlog_out(DEBUG4) << "Difference in +P = " << dPplus << "  -P = " << dPminus << endmsg;

      TrackExtended * combinedTrack = CombineTracks(firstTrackExt,secondTrackExt, _maxAllowedPercentageOfOutliersForTrackCombination, true);
      
      if(combinedTrack != NULL){
        streamlog_out(DEBUG4) << "CombinedTrack " << combinedTrack->getNDF() << " c.f. " << firstTrackExt->getNDF()+secondTrackExt->getNDF()+5 << endmsg;
      }else{
        streamlog_out(DEBUG4) << "Could not combine track " << endmsg;
      }
      delete combinedTrack->getGroupTracks();
      delete combinedTrack;
      streamlog_out(DEBUG4) << " Overlap = " << SegmentRadialOverlap(firstTrackExt, secondTrackExt) << " veto = " << VetoMerge(firstTrackExt, secondTrackExt) << endmsg;
      
      streamlog_out(DEBUG4) << endmsg;
      
    }
  }
  
  catch( ... ){};
  
}    
*/

/*
void FullLDCTrackingAlg::GeneralSorting(int * index, float * val, int direct, int nVal) {
   // Sorting of index vector in ascending (0) /descending (!=0) order of val
  
  float valOne, valTwo, valTemp;
  int   indTemp;
  for (int i=0; i<nVal; ++i) {
    index[i] = i;
  }
  
  for (int i = 0 ; i < nVal-1; i++) {
    for (int j = 0; j < nVal-i-1; j++) {      
      valOne = val[j];
      valTwo = val[j+1];
      bool order = valOne > valTwo;
      if (direct>0) 
        order = valOne <= valTwo;
      if( order )
      {
        valTemp = val[j];
        val[j] = val[j+1];
        val[j+1] = valTemp;
        indTemp = index[j];
        index[j] = index[j+1];
        index[j+1] = indTemp;
      }
    }  
  }
  
  
}
*/

int FullLDCTrackingAlg::SegmentRadialOverlap(TrackExtended* first, TrackExtended* second){
  
  
  int nTrkGrpFirst = 0;
  int nTrkGrpSecond = 0;
  ConstTrackerHitVec hitvecFirst;
  ConstTrackerHitVec hitvecSecond;
  GroupTracks * groupFirst = first->getGroupTracks();
  GroupTracks * groupSecond = second->getGroupTracks();
  
  if(groupFirst!=NULL){
    
    TrackExtendedVec tracksInGroupFirst = groupFirst->getTrackExtendedVec();
    nTrkGrpFirst = int(tracksInGroupFirst.size());
    
    for (int iTrkGrp=0;iTrkGrp<nTrkGrpFirst;++iTrkGrp) {
      
      TrackExtended * trkGrp = tracksInGroupFirst[iTrkGrp];
      TrackerHitExtendedVec hitVec = trkGrp->getTrackerHitExtendedVec();
      
      for(unsigned int i =0; i<hitVec.size(); ++i){
        hitvecFirst.push_back(hitVec[i]->getTrackerHit());        
      }
    }
  }
  if(groupSecond!=NULL){
    
    TrackExtendedVec tracksInGroupSecond = groupSecond->getTrackExtendedVec();
    nTrkGrpSecond = int(tracksInGroupSecond.size());
    
    for (int iTrkGrp=0;iTrkGrp<nTrkGrpSecond;++iTrkGrp) {
      TrackExtended * trkGrp = tracksInGroupSecond[iTrkGrp];
      TrackerHitExtendedVec hitVec = 
      trkGrp->getTrackerHitExtendedVec();
      
      for(unsigned int i=0;i<hitVec.size();++i){
        hitvecSecond.push_back(hitVec[i]->getTrackerHit());
      }
    }
  }
  
  
  int nhitsFirst = (int)hitvecFirst.size();
  int nhitsSecond = (int)hitvecSecond.size();
  int count = 0;
  for(int i =0;i<nhitsFirst;++i){
    float xi = (float)hitvecFirst[i].getPosition()[0];
    float yi = (float)hitvecFirst[i].getPosition()[1];
    float ri = sqrt(xi*xi+yi*yi);
    if(ri < _tpc_inner_r || ri > _tpc_pad_height)continue;
    for(int j =0;j<nhitsSecond;++j){
      float xj = (float)hitvecSecond[j].getPosition()[0];
      float yj = (float)hitvecSecond[j].getPosition()[1];
      float rj = sqrt(xj*xj+yj*yj);
      if(fabs(ri-rj)<_tpc_pad_height/2.0)count++;
    }
  }  
  return count;
}

/*
 
 veto merger if the momentum of either track is less than 2.5 GeV
 or if following a full fit the NDF+10 of the combined tracks is less than the NDF_first + NDF_second
 
 NOTE: This function will try a full fit using CombineTracks if the momentum of both tracks is greater than VetoMergeMomentumCut
 
 */

bool FullLDCTrackingAlg::VetoMerge(TrackExtended* firstTrackExt, TrackExtended* secondTrackExt){
    
  debug() << "FullLDCTrackingAlg::VetoMerge called for " << firstTrackExt << " and " << secondTrackExt << endmsg;
  
  const float d0First = firstTrackExt->getD0();
  const float z0First = firstTrackExt->getZ0();
  const float omegaFirst = firstTrackExt->getOmega();
  const float tanLFirst = firstTrackExt->getTanLambda();
  const float phi0First = firstTrackExt->getPhi();
  
  const float d0Second = secondTrackExt->getD0();
  const float z0Second = secondTrackExt->getZ0();
  const float omegaSecond = secondTrackExt->getOmega();
  const float tanLSecond = secondTrackExt->getTanLambda();
  const float phi0Second = secondTrackExt->getPhi();        
  
  HelixClass firstHelix;
  firstHelix.Initialize_Canonical(phi0First,d0First,z0First,omegaFirst,tanLFirst,_bField);
  const float pxFirst = firstHelix.getMomentum()[0];
  const float pyFirst = firstHelix.getMomentum()[1];
  const float pzFirst = firstHelix.getMomentum()[2];        
  
  HelixClass secondHelix;
  secondHelix.Initialize_Canonical(phi0Second,d0Second,z0Second,omegaSecond,tanLSecond,_bField);
  const float pxSecond = secondHelix.getMomentum()[0];
  const float pySecond = secondHelix.getMomentum()[1];
  const float pzSecond = secondHelix.getMomentum()[2];      
  const float pFirst = sqrt(pxFirst*pxFirst+pyFirst*pyFirst+pzFirst*pzFirst);
  const float pSecond = sqrt(pxSecond*pxSecond+pySecond*pySecond+pzSecond*pzSecond);
  
  if(pFirst<_vetoMergeMomentumCut || pSecond<_vetoMergeMomentumCut) {
    debug() << "FullLDCTrackingAlg::VetoMerge do not veto as below momentum cut of 2.5 : pFirst = " << pFirst << " pSecond = " << pSecond << endmsg;
    return false;
  }

  bool veto = false;
  
  bool testCombinationOnly=true;
  TrackExtended * combinedTrack = CombineTracks(firstTrackExt,secondTrackExt,_maxAllowedPercentageOfOutliersForTrackCombination,testCombinationOnly);
  
  if(combinedTrack!=NULL){
    //SJA:FIXME hardcoded cut: here the check is that no more than 7 hits have been rejected in the combined fit.
    if( combinedTrack->getNDF()+15 < firstTrackExt->getNDF() + secondTrackExt->getNDF()+5 ) {
      debug() << "FullLDCTrackingAlg::VetoMerge fails NDF cut " << endmsg;
      veto=true ;
    }
  
    delete combinedTrack->getGroupTracks();
    delete combinedTrack;

  }
  else {
    debug() << "FullLDCTrackingAlg::VetoMerge fails CombineTracks(firstTrackExt,secondTrackExt,true) test" << endmsg;
    veto = true;
  }
  
  if(SegmentRadialOverlap(firstTrackExt,secondTrackExt)>10) {
    debug() << "FullLDCTrackingAlg::VetoMerge fails SegmentRadialOverlap test " << endmsg;
    veto=true;
  }

  return veto;
}

StatusCode FullLDCTrackingAlg::finalize() { 
  delete _encoder ;
  return StatusCode::SUCCESS;
}

void FullLDCTrackingAlg::setupGearGeom(){
  
  auto _gear = service<IGearSvc>("GearSvc");
  gearMgr = _gear->getGearMgr();

  _bField = gearMgr->getBField().at( gear::Vector3D( 0.,0.,0.)  ).z() ;
  
  //-- VXD Parameters--
  _nLayersVTX = 0 ;
  const gear::VXDParameters* pVXDDetMain = 0;
  const gear::VXDLayerLayout* pVXDLayerLayout = 0;
  
  try{
    debug() << " filling VXD parameters from gear::SITParameters " << endmsg;
    
    pVXDDetMain = &(gearMgr->getVXDParameters());
    pVXDLayerLayout = &(pVXDDetMain->getVXDLayerLayout());
    _nLayersVTX = pVXDLayerLayout->getNLayers();
  }
  catch( ... ){
    debug() << " ### gear::VXDParameters Not Present in GEAR FILE" << endmsg;
  }
    
  //-- SIT Parameters--
  _nLayersSIT = 0 ;
  const gear::ZPlanarParameters* pSITDetMain = 0;
  const gear::ZPlanarLayerLayout* pSITLayerLayout = 0;
  
  try{
    debug() << " filling SIT parameters from gear::SITParameters " << endmsg;
    
    pSITDetMain = &(gearMgr->getSITParameters());
    pSITLayerLayout = &(pSITDetMain->getZPlanarLayerLayout());
    _nLayersSIT = pSITLayerLayout->getNLayers();
  }
  catch( ... ){
    debug() << " ### gear::SITParameters Not Present in GEAR FILE" << endmsg;
  }
  
  if( _nLayersSIT == 0 ){
    // try the old LOI style key value pairs as defined in the SSit03 Mokka drive
    try{
      debug() << "  FullLDCTrackingAlg - Simple Cylinder Based SIT using parameters defined by SSit03 Mokka driver " << endmsg;
      
      // SIT
      const gear::GearParameters& pSIT = gearMgr->getGearParameters("SIT");
      
      const EVENT::DoubleVec& SIT_r   =  pSIT.getDoubleVals("SITLayerRadius" )  ;
      const EVENT::DoubleVec& SIT_hl  =  pSIT.getDoubleVals("SITSupportLayerHalfLength" )  ;
      
      _nLayersSIT = SIT_r.size() ; 
      
      if (_nLayersSIT != SIT_r.size() || _nLayersSIT != SIT_hl.size()) {
	fatal() << "FullLDCTrackingAlg Miss-match between DoubleVec and nlayers exit(1) called from file " << __FILE__ << " line " << __LINE__  << endmsg;
        exit(1);
      }
    }
    catch( ... ){
      debug() << " ### gear::SIT Parameters from as defined in SSit03 Not Present in GEAR FILE" << endmsg;
    } 
  }
  
  //-- SET Parameters--
  _nLayersSET = 0 ;
  const gear::ZPlanarParameters* pSETDetMain = 0;
  const gear::ZPlanarLayerLayout* pSETLayerLayout = 0;
  
  try{
    debug() << " filling SET parameters from gear::SETParameters " << endmsg;
    
    pSETDetMain = &(gearMgr->getSETParameters());
    pSETLayerLayout = &(pSETDetMain->getZPlanarLayerLayout());
    _nLayersSET = pSETLayerLayout->getNLayers();
  }
  catch( ... ){
    debug() << " ### gear::SETParameters Not Present in GEAR FILE" << endmsg;
  }
  
  if( _nLayersSET == 0 ){
    // try the old LOI style key value pairs as defined in the SSet02 Mokka drive
    try{
      debug() << "  FullLDCTrackingAlg - Simple Cylinder Based SET using parameters defined by SSet02 Mokka driver " << endmsg;
      
      // SET
      const gear::GearParameters& pSET = gearMgr->getGearParameters("SET");
      
      const EVENT::DoubleVec& SET_r   =  pSET.getDoubleVals("SETLayerRadius" )  ;
      const EVENT::DoubleVec& SET_hl  =  pSET.getDoubleVals("SETSupportLayerHalfLength" )  ;
      
      _nLayersSET = SET_r.size() ;
      
      if (_nLayersSET != SET_r.size() || _nLayersSET != SET_hl.size()) {
	fatal() << "FullLDCTrackingAlg Miss-match between DoubleVec and nlayers exit(1) called from file " << __FILE__ << " line " << __LINE__  << endmsg;
        exit(1);
      }
    }
    catch( ... ){
      debug() << " ### gear::SET Parameters from as defined in SSet02 Not Present in GEAR FILE" << endmsg;
    }
  }
    
  //-- FTD Parameters--
  _petalBasedFTDWithOverlaps = false;  
  _nLayersFTD = 0;
  
  try{
    debug() << " filling FTD parameters from gear::FTDParameters " << endmsg;
    
    const gear::FTDParameters&   pFTD      =gearMgr->getFTDParameters();
    const gear::FTDLayerLayout&  ftdlayers = pFTD.getFTDLayerLayout() ;
    
    _nLayersFTD = ftdlayers.getNLayers() ;
    
    for (unsigned int disk=0; disk < _nLayersFTD; ++disk) {
      
      _zLayerFTD.push_back( ftdlayers.getSensitiveZposition(disk, 0, 1) ); // front petal even numbered
      
      if ( ftdlayers.getNPetals(disk) > 0) {
        _zLayerFTD.push_back( ftdlayers.getSensitiveZposition(disk, 1, 1) );  // front petal odd numbered
        _petalBasedFTDWithOverlaps = true;
      }
      
    }
    
    // SJA: Here we increase the size of _nlayersFTD as we are treating the 
    _nLayersFTD =_zLayerFTD.size() ;     
    
  }
  catch( ... ){
    debug() << " ### gear::FTDParameters Not Present in GEAR FILE" << endmsg;
  } 
  
  if( _nLayersFTD == 0 ){
    // FTD
    try{
      debug() << "  FullLDCTrackingAlg - Simple Disc Based FTD using parameters defined by SFtd05 Mokka driver " << endmsg;
      
      const gear::GearParameters& pFTD = gearMgr->getGearParameters("FTD");
      
      const EVENT::DoubleVec* pFTD_z   = NULL;
      
      pFTD_z = &pFTD.getDoubleVals("FTDZCoordinate" )  ;
      
      _nLayersFTD = pFTD_z->size();
      
      for (unsigned int i = 0; i<_nLayersFTD; ++i) {
        _zLayerFTD.push_back((*pFTD_z)[i]);
      }
    }
    catch( ... ){
      debug() << " ### gear::FTD Parameters as defined in SFtd05 Not Present in GEAR FILE" << endmsg;
    } 
  }
}
