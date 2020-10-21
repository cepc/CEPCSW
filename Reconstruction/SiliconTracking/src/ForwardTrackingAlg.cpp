#include "ForwardTrackingAlg.h"
#include "GearSvc/IGearSvc.h"
#include "TrackSystemSvc/ITrackSystemSvc.h"
#include "DataHelper/Navigation.h"

#include "edm4hep/TrackerHit.h"
#include "edm4hep/TrackerHitConst.h"
#include "edm4hep/Track.h"

#include "UTIL/ILDConf.h"

//#include "MarlinCED.h"

#include "gear/GEAR.h"
#include "gear/GearParameters.h"
#include "gear/BField.h"
#include "gear/FTDParameters.h"
#include "gear/FTDLayerLayout.h"

//----From KiTrack-----------------------------
#include "KiTrack/SubsetHopfieldNN.h"
#include "KiTrack/SubsetSimple.h"
#include "KiTrack/SegmentBuilder.h"
#include "KiTrack/Automaton.h"

//----From KiTrackMarlin-----------------------
#include "ILDImpl/FTDTrack.h"
#include "ILDImpl/FTDHit01.h"
#include "ILDImpl/FTDNeighborPetalSecCon.h"
#include "ILDImpl/FTDSectorConnector.h"
#include "Tools/KiTrackMarlinTools.h"
//#include "Tools/KiTrackMarlinCEDTools.h"
#include "Tools/FTDHelixFitter.h"

using namespace MarlinTrk ;

// Used to fedine the quality of the track output collection
const int ForwardTrackingAlg::_output_track_col_quality_GOOD = 1;
const int ForwardTrackingAlg::_output_track_col_quality_FAIR = 2;
const int ForwardTrackingAlg::_output_track_col_quality_POOR = 3;

DECLARE_COMPONENT( ForwardTrackingAlg )

ForwardTrackingAlg::ForwardTrackingAlg(const std::string& name, ISvcLocator* svcLoc)
: GaudiAlgorithm(name, svcLoc) {
  //_description = "ForwardTracking reconstructs tracks through the FTD" ;

  declareProperty("FTDPixelHitCollection", _inFTDPixelColHdl, "Handle of the Input FTD TrackerHits collection");
  declareProperty("FTDSpacePointCollection", _inFTDSpacePointColHdl, "Handle of the Input FTD SpacePoints collection");
  declareProperty("FTDRawHitCollection", _inFTDRawColHdl, "Handle of the FTD SpacePoints raw hit collection collection");
  //std::vector< std::string > collections;
  //collections.push_back( "FTDTrackerHits" );
  //collections.push_back( "FTDSpacePoints" );
  declareProperty("ForwardTrackCollection", _outColHdl, "Handle of the ForwarTrack output collection");
}

StatusCode ForwardTrackingAlg::initialize(){
  debug() << "   init called  " << endmsg;
  
  _nRun = 0 ;
  _nEvt = 0 ;

  _useCED = false; // Setting this to on will initialise CED in the processor and tracks or segments (from the CA)
  // can be printed. As this is mainly used for debugging it is not a steerable parameter.
  //if( _useCED )MarlinCED::init(this) ;    //CED

  // Now set min and max values for all the criteria
  for( unsigned i=0; i < _criteriaNames.size(); i++ ){
    std::vector< float > emptyVec;
    _critMinima[_criteriaNames[i]] = emptyVec;
    _critMaxima[_criteriaNames[i]] = emptyVec;
  }
  if(_critMinimaInit.size()!=_critMaximaInit.size()){
    warning() << "number of Criteria min values != max values, will be filled as 0 for less" << endmsg;
  }
  unsigned nCritValues = std::max(_critMinimaInit.size(),_critMaximaInit.size());
  //if(nCritValues>_criteriaNames.size()){
  //  warning() << "number of Criteria values > name's, will be discarded" << endmsg;
  //}
  //nCritValues = std::min(nCritValues, _criteriaNames.size());
  for(unsigned i=0; i<nCritValues; i++){
    int iCritName = i%_criteriaNames.size();
    std::string critName = _criteriaNames[iCritName];
    if(i<_critMinimaInit.size()) _critMinima[critName].push_back(_critMinimaInit[i]);
    else                         _critMinima[critName].push_back(0.);
    if(i<_critMaximaInit.size()) _critMaxima[critName].push_back(_critMaximaInit[i]);
    else                         _critMaxima[critName].push_back(0.);
  }
  debug() << "Criteria table:" << endmsg;
  for(unsigned i=0; i<_criteriaNames.size(); i++ ){
    std::string critName = _criteriaNames[i];
    debug() << "---- " << critName;
    for(unsigned j=0; j<_critMinima[critName].size(); j++){
      debug() << " [" << _critMinima[critName][j] << ", " << _critMaxima[critName][j] << "]"; 
    }
    debug() << endmsg;
  }

  auto _gear = service<IGearSvc>("GearSvc");
  if ( !_gear ) {
    error() << "Failed to find GearSvc ..." << endmsg;
    return StatusCode::FAILURE;
  }
  gear::GearMgr* gearMgr = _gear->getGearMgr();
     
  /**********************************************************************************************/
  /*       Make a SectorSystemFTD                                                               */
  /**********************************************************************************************/

  // The SectorSystemFTD is the object translating the sectors of the hits into layers, modules etc. and vice versa
  const gear::FTDParameters& ftdParams = gearMgr->getFTDParameters() ;
  const gear::FTDLayerLayout& ftdLayers = ftdParams.getFTDLayerLayout() ;
  int nLayers = ftdLayers.getNLayers() + 1; // we add one layer for the IP
  int nModules = ftdLayers.getNPetals(0);
  int nSensors = ftdLayers.getNSensors(0);
  
  // make sure we take the highest number of modules / sensors available
  for( int i=1; i < nLayers - 1; i++){
     
    if( ftdLayers.getNPetals(i) > nModules ) nModules = ftdLayers.getNPetals(i); 
    if( ftdLayers.getNSensors(i) > nSensors ) nSensors = ftdLayers.getNSensors(i);
    
  }
  
  debug() << "SectorSystemFTD is using " << nLayers << " layers (including one for the IP), " << nModules << " petals and " << nSensors << " sensors." << endmsg;
   
  _sectorSystemFTD = new SectorSystemFTD( nLayers, nModules , nSensors );

  // Get the B Field in z direction
  _Bz = gearMgr->getBField().at( gear::Vector3D(0., 0., 0.) ).z();    //The B field in z direction
  
  /**********************************************************************************************/
  /*       Initialise the MarlinTrkSystem, needed by the tracks for fitting                     */
  /**********************************************************************************************/
  auto _trackSystemSvc = service<ITrackSystemSvc>("TrackSystemSvc");
  if ( !_trackSystemSvc ) {
    error() << "Failed to find TrackSystemSvc ..." << endmsg;
    return StatusCode::FAILURE;
  }
  _trkSystem =  _trackSystemSvc->getTrackSystem(this);

  if( _trkSystem == 0 ){
    error() << "Cannot initialize MarlinTrkSystem of Type: KalTest" <<endmsg;
    return StatusCode::FAILURE;
  }

  // set the options   
  _trkSystem->setOption( MarlinTrk::IMarlinTrkSystem::CFG::useQMS,        _MSOn ) ;       //multiple scattering
  _trkSystem->setOption( MarlinTrk::IMarlinTrkSystem::CFG::usedEdx,       _ElossOn) ;     //energy loss
  _trkSystem->setOption( MarlinTrk::IMarlinTrkSystem::CFG::useSmoothing,  _SmoothOn) ;    //smoothing

  // initialise the tracking system
  _trkSystem->init() ;

  /**********************************************************************************************/
  /*       Do a few checks, if the set parameters are right                                     */
  /**********************************************************************************************/
    
  // Only use allowed methods to find subsets. 
  assert( ( _bestSubsetFinder == "None" ) || ( _bestSubsetFinder == "SubsetHopfieldNN" ) || ( _bestSubsetFinder == "SubsetSimple" ) );
  
  // Use a sensible chi2prob cut. (chi squared probability, like any probability must range from 0 to 1)
  assert( _chi2ProbCut >= 0. );
  assert( _chi2ProbCut <= 1. );
     
  // Make sure, every used criterion exists and has at least one min and max set
  for( unsigned i=0; i<_criteriaNames.size(); i++ ){
    
    std::string critName = _criteriaNames[i];
    
    ICriterion* crit = Criteria::createCriterion( critName ); //throws an exception if the criterion is non existent
    delete crit;
    
    assert( !_critMinima[ critName ].empty() );
    assert( !_critMaxima[ critName ].empty() );
    
  }
  return GaudiAlgorithm::initialize();
}

StatusCode ForwardTrackingAlg::execute(){
  debug() << " processing event number " << _nEvt << endmsg;

  auto trkCol = _outColHdl.createAndPut();
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //                                                                                                              //
  //                                 ForwardTracking                                                              //
  //                                                                                                              //
  //                            Track Reconstruction in the FTD                                                   //
  //                                                                                                              //
  //                                                                                                              //
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   
  //--CED (only used for debugging )---------------------------------------
  // Reset drawing buffer and START drawing collection
  /*
  if( _useCED ){
    MarlinCED::newEvent(this , 0) ; 
    CEDPickingHandler &pHandler=CEDPickingHandler::getInstance();
    pHandler.update(evt); 
  }
  */
  //-----------------------------------------------------------------------
    
  // Reset the quality flag of the output track collection (we start with the assumption that our results are good.
  // If anything happens along the way, we modify this value )
  _output_track_col_quality = _output_track_col_quality_GOOD;
  
  std::vector< IHit* > hitsTBD; //Hits to be deleted at the end
  _map_sector_hits.clear();
     
  /**********************************************************************************************/
  /*    Read in the collections, create hits from the TrackerHits and store them in a map       */
  /**********************************************************************************************/
  
  debug() << "\t\t---Reading in Collections---" << endmsg;
  Navigation::Instance()->Initialize();
  std::vector<const edm4hep::TrackerHitCollection*> hitFTDCollections;
  int pixelCollectionID = -1; 
  try {
    auto hitFTDPixelCol = _inFTDPixelColHdl.get();
    pixelCollectionID = hitFTDPixelCol->getID();
    hitFTDCollections.push_back(hitFTDPixelCol);
    Navigation::Instance()->AddTrackerHitCollection(hitFTDPixelCol);
  }
  catch ( GaudiException &e ) {
    debug() << "Collection " << _inFTDPixelColHdl.fullKey() << " is unavailable in event " << _nEvt << endmsg;
  }

  try {
    auto hitFTDSpacePointCol = _inFTDSpacePointColHdl.get();
    hitFTDCollections.push_back(hitFTDSpacePointCol);
    Navigation::Instance()->AddTrackerHitCollection(hitFTDSpacePointCol);
    //const edm4hep::TrackerHitCollection* rawHitCol = nullptr;
    try{
      auto rawHitCol = _inFTDRawColHdl.get();
      Navigation::Instance()->AddTrackerHitCollection(rawHitCol);
    }
    catch ( GaudiException &e ) {
      fatal() << "Collection " << _inFTDRawColHdl.fullKey() << " is unavailable in event " << _nEvt << endmsg;
    }
  }
  catch ( GaudiException &e ) {
    debug() << "Collection " << _inFTDSpacePointColHdl.fullKey() << " is unavailable in event " << _nEvt << endmsg;
  }

  if(hitFTDCollections.size()==0){
    return StatusCode::SUCCESS;
  }
  /*
  const edm4hep::TrackerHitCollection* rawHitCol = nullptr;
  if(1){
    try{
      rawHitCol = _inFTDRawColHdl.get();
    }
    catch ( GaudiException &e ) {
      fatal() << "Collection " << _inFTDRawColHdl.fullKey() << " is unavailable in event " << _nEvt << endmsg;
    }
  }
  */
  for( unsigned iCol=0; iCol < hitFTDCollections.size(); iCol++ ){ //read in all input collections
    unsigned nHits = hitFTDCollections[iCol]->size();
    debug() << "Number of hits in collection " << hitFTDCollections[iCol]->getID() << ": " << nHits << endmsg;

    for(auto trackerHit : *hitFTDCollections[iCol]){
      if(pixelCollectionID==hitFTDCollections[iCol]->getID()){
	if ( UTIL::BitSet32( trackerHit.getType() )[ UTIL::ILDTrkHitTypeBit::ONE_DIMENSIONAL ] ) continue;
      }
      edm4hep::ConstTrackerHit hit = trackerHit;
      debug() << "hit " << trackerHit.id() << " " << KiTrackMarlin::getCellID0Info( trackerHit.getCellID() ) 
	      << " " << KiTrackMarlin::getPositionInfo( hit )<< endmsg;
         
      //Make an FTDHit01 from the TrackerHit 
      FTDHit01* ftdHit = new FTDHit01 ( trackerHit , _sectorSystemFTD );
      hitsTBD.push_back(ftdHit); //so we can easily delete every created hit afterwards
      
      _map_sector_hits[ ftdHit->getSector() ].push_back( ftdHit );         
    }
  }
     
  if( !_map_sector_hits.empty() ){
    /**********************************************************************************************/
    /*                Check if no sector is overflowing with hits                                 */
    /**********************************************************************************************/
    std::map< int , std::vector< IHit* > >::iterator it;
      
    for( it=_map_sector_hits.begin(); it != _map_sector_hits.end(); it++ ){
      int nHits = it->second.size();
      debug() << "Number of hits in sector " << it->first << " = " << nHits << endmsg;
      
      if( nHits > _maxHitsPerSector ){
	
	it->second.clear(); //delete the hits in this sector, it will be dropped
        
	error() << " ### EVENT " << _nEvt << " :: RUN " << _nRun << " \n ### Number of Hits in FTD Sector " << it->first << ": " << nHits << " > "
		<< _maxHitsPerSector << " (MaxHitsPerSector)\n : This sector will be dropped from track search, and QualityCode set to \"Poor\" " << endmsg;
	
	_output_track_col_quality = _output_track_col_quality_POOR; // We had to drop hits, so the quality of the result is decreased
      }
    }
      
    /**********************************************************************************************/
    /*                Check the possible connections of hits on overlapping petals                */
    /**********************************************************************************************/
    debug() << "\t\t---Overlapping Hits---" << endmsg;
      
    std::map< IHit* , std::vector< IHit* > > map_hitFront_hitsBack = getOverlapConnectionMap( _map_sector_hits, _sectorSystemFTD, _overlappingHitsDistMax);

    /**********************************************************************************************/
    /*                Add the IP as virtual hit for forward and backward                          */
    /**********************************************************************************************/
    IHit* virtualIPHitForward = createVirtualIPHit(1 , _sectorSystemFTD );
    hitsTBD.push_back( virtualIPHitForward );
    _map_sector_hits[ virtualIPHitForward->getSector() ].push_back( virtualIPHitForward );
    
    IHit* virtualIPHitBackward = createVirtualIPHit(-1 , _sectorSystemFTD );
    hitsTBD.push_back( virtualIPHitBackward );
    _map_sector_hits[ virtualIPHitBackward->getSector() ].push_back( virtualIPHitBackward );
     
    /**********************************************************************************************/
    /*                SegmentBuilder and Cellular Automaton                                       */
    /**********************************************************************************************/
    unsigned round = 0; // the round we are in
    std::vector < RawTrack > rawTracks;

    // The following while loop ideally only runs once. (So we do round 0 and everything works)
    // It will repeat as long as the Automaton creates too many connections and as long as there are new criteria
    // parameters to use to cut down the problem.
    // Ideally already in round 0, there is a reasonable number of connections (not more than _maxConnectionsAutomaton), 
    // so the loop will be left. If however there are too many connections we stay in the loop and use 
    // (hopefully) tighter cut offs (if provided in the steering). This should prevent combinatorial breakdown
    // for very evil events.
    while( setCriteria( round ) ){
      round++; // count up the round we are in
               
      /**********************************************************************************************/
      /*                Build the segments                                                          */
      /**********************************************************************************************/
      debug() << "\t\t---SegementBuilder---" << endmsg;
         
      //Create a segmentbuilder
      SegmentBuilder segBuilder( _map_sector_hits );
      
      segBuilder.addCriteria ( _crit2Vec ); // Add the criteria on when to connect two hits. The vector has been filled by the method setCriteria
         
      //Also load hit connectors
      unsigned layerStepMax = 1; // how many layers to go at max
      unsigned petalStepMax = 1; // how many petals to go at max
      unsigned lastLayerToIP = 5;// layer 1,2,3 and 4 get connected directly to the IP
      FTDSectorConnector secCon( _sectorSystemFTD , layerStepMax , petalStepMax , lastLayerToIP );
               
      segBuilder.addSectorConnector ( & secCon ); // Add the sector connector (so the SegmentBuilder knows what hits from different sectors it is allowed to look for connections)
               
      // And get out the Cellular Automaton with the 1-segments 
      Automaton automaton = segBuilder.get1SegAutomaton();
      
      // Check if there are not too many connections
      if( automaton.getNumberOfConnections() > unsigned( _maxConnectionsAutomaton ) ){
	debug() << "Redo the Automaton with different parameters, because there are too many connections:" << endmsg
		<< "\tconnections( " << automaton.getNumberOfConnections() << " ) > MaxConnectionsAutomaton( " << _maxConnectionsAutomaton << " )" << endmsg;
	continue;
      }
               
      /**********************************************************************************************/
      /*                Automaton                                                                   */
      /**********************************************************************************************/
      debug() << "\t\t---Automaton---" << endmsg;
      
      //if( _useCED ) KiTrackMarlin::drawAutomatonSegments( automaton ); // draws the 1-segments (i.e. hits)
               
      /*******************************/
      /*      2-hit segments         */
      /*******************************/
      debug() << "\t\t--2-hit-Segments--" << endmsg;
         
      //debug() << "Automaton has " << automaton.getTracks( 3 ).size() << " track candidates" << endmsg; //should be commented out, because it takes time
         
      automaton.clearCriteria();
      automaton.addCriteria( _crit3Vec );  // Add the criteria for 3 hits (i.e. 2 2-hit segments )
               
      // Let the automaton lengthen its 1-hit-segments to 2-hit-segments
      automaton.lengthenSegments();
               
      // So now we have 2-hit-segments and are ready to perform the Cellular Automaton.
      
      // Perform the automaton
      automaton.doAutomaton();
               
      // Clean segments with bad states
      automaton.cleanBadStates();
              
      // Reset the states of all segments
      automaton.resetStates();
      
      //debug() << "Automaton has " << automaton.getTracks( 3 ).size() << " track candidates" << endmsg; //should be commented out, because it takes time
                  
      // Check if there are not too many connections
      if( automaton.getNumberOfConnections() > unsigned( _maxConnectionsAutomaton ) ){
	debug() << "Redo the Automaton with different parameters, because there are too many connections:" << endmsg
		<< "\tconnections( " << automaton.getNumberOfConnections() << " ) > MaxConnectionsAutomaton( " << _maxConnectionsAutomaton << " )" << endmsg;
	continue;
      }
         
      /*******************************/
      /*      3-hit segments         */
      /*******************************/
      debug() << "\t\t--3-hit-Segments--" << endmsg;
         
      automaton.clearCriteria();
      automaton.addCriteria( _crit4Vec );      
      
      // Lengthen the 2-hit-segments to 3-hits-segments
      automaton.lengthenSegments();
      
      // Perform the Cellular Automaton
      automaton.doAutomaton();
      
      //Clean segments with bad states
      automaton.cleanBadStates();
               
      //Reset the states of all segments
      automaton.resetStates();
                  
      //debug() << "Automaton has " << automaton.getTracks( 3 ).size() << " track candidates" << endmsg; //should be commented out, because it takes time
      
      // Check if there are not too many connections
      if( automaton.getNumberOfConnections() > unsigned( _maxConnectionsAutomaton ) ){
	debug() << "Redo the Automaton with different parameters, because there are too many connections:" << endmsg
		<< "\tconnections( " << automaton.getNumberOfConnections() << " ) > MaxConnectionsAutomaton( " << _maxConnectionsAutomaton << " )" << endmsg;
	continue;
      }
         
      // get the raw tracks (raw track = just a vector of hits, the most rudimentary form of a track)
      rawTracks = automaton.getTracks( 3 );
      
      break; // if we reached this place all went well and we don't need another round --> exit the loop
    }
      
    debug() << "Automaton returned " << rawTracks.size() << " raw tracks " << endmsg;
        
    /**********************************************************************************************/
    /*                Add the overlapping hits                                                    */
    /**********************************************************************************************/
    debug() << "\t\t---Add hits from overlapping petals + fit + helix and Kalman cuts---" << endmsg;
    
    std::vector <ITrack*> trackCandidates;
    
    // for all raw tracks we got from the automaton
    for( unsigned i=0; i < rawTracks.size(); i++){
      RawTrack rawTrack = rawTracks[i];
      
      _nTrackCandidates++;
              
      // get all versions of the track plus hits from overlapping petals
      std::vector < RawTrack > rawTracksPlus = getRawTracksPlusOverlappingHits( rawTrack, map_hitFront_hitsBack );
      
      debug() << "For raw track number " << i << " there are " << rawTracksPlus.size() << " versions" << endmsg;
                  
      /**********************************************************************************************/
      /*                Make track candidates, fit them and throw away bad ones                     */
      /**********************************************************************************************/
      std::vector< ITrack* > overlappingTrackCands;
        
      for( unsigned j=0; j < rawTracksPlus.size(); j++ ){
	_nTrackCandidatesPlus++;
	    
	RawTrack rawTrackPlus = rawTracksPlus[j];
        
	if( rawTrackPlus.size() < unsigned( _hitsPerTrackMin ) ){
	  debug() << "Trackversion discarded, too few hits: only " << rawTrackPlus.size() << " < " << _hitsPerTrackMin << "(hitsPerTrackMin)" << endmsg;
	  continue;
	}
            
	FTDTrack* trackCand = new FTDTrack( _trkSystem );
            
	// add the hits to the track
	for( unsigned k=0; k<rawTrackPlus.size(); k++ ){
	  IFTDHit* ftdHit = dynamic_cast< IFTDHit* >( rawTrackPlus[k] ); // cast to IFTDHits, as needed for an FTDTrack
	  if( ftdHit != NULL ) trackCand->addHit( ftdHit );
	  else debug() << "Hit " << rawTrackPlus[k] << " could not be casted to IFTDHit" << endmsg;
	}
            
	std::vector< IHit* > trackCandHits = trackCand->getHits();
	debug() << "Fitting track candidate with " << trackCandHits.size() << " hits" << endmsg;
        
	for( unsigned k=0; k < trackCandHits.size(); k++ ) debug() << trackCandHits[k]->getPositionInfo();
	debug() << endmsg;
            
	/*-----------------------------------------------*/
	/*                Helix Fit                      */
	/*-----------------------------------------------*/
	debug() << "Fitting with Helix Fit" << endmsg;
	try{
	  FTDHelixFitter helixFitter( trackCand->getLcioTrack() );
	  float chi2OverNdf = helixFitter.getChi2() / float( helixFitter.getNdf() );
	  debug() << "chi2OverNdf = " << chi2OverNdf << endmsg;
          
	  if( chi2OverNdf > _helixFitMax ){
	    debug() << "Discarding track because of bad helix fit: chi2/ndf = " << chi2OverNdf << endmsg;
	    delete trackCand;
	    continue;
	  }
	  else debug() << "Keeping track because of good helix fit: chi2/ndf = " << chi2OverNdf << endmsg;
	}
	catch( FTDHelixFitterException& e ){
	  debug() << "Track rejected, because fit failed: " <<  e.what() << endmsg;
	  delete trackCand;
	  continue;
	}
            
	/*-----------------------------------------------*/
	/*                Kalman Fit                      */
	/*-----------------------------------------------*/
        
	debug() << "Fitting with Kalman Filter" << endmsg;
	try{
	  trackCand->fit();
                  
	  debug() << " Track " << trackCand 
		  << " chi2Prob = " << trackCand->getChi2Prob() 
		  << "( chi2=" << trackCand->getChi2() 
		  <<", Ndf=" << trackCand->getNdf() << " )" << endmsg;
	                    
	  if ( trackCand->getChi2Prob() >= _chi2ProbCut ){
	    debug() << "Track accepted (chi2prob " << trackCand->getChi2Prob() << " >= " << _chi2ProbCut << endmsg;
	  }
	  else{
	    debug() << "Track rejected (chi2prob " << trackCand->getChi2Prob() << " < " << _chi2ProbCut << endmsg;
	    delete trackCand;
            
	    continue;
	  }
	}
	catch( FitterException& e ){
	  debug() << "Track rejected, because fit failed: " <<  e.what() << endmsg;
	  delete trackCand;
	  continue;
	}
            
	// If we reach this point than the track got accepted by all cuts
	overlappingTrackCands.push_back( trackCand );
      }
         
      /**********************************************************************************************/
      /*                Take the best version of the track                                          */
      /**********************************************************************************************/
      // Now we have all versions of one track, coming from adding possible hits from overlapping petals.
      if( _takeBestVersionOfTrack ){ // we want to take only the best version
	debug() << "Take the version of the track with best quality from " << overlappingTrackCands.size() << " track candidates" << endmsg;
            
	if( !overlappingTrackCands.empty() ){
	  ITrack* bestTrack = overlappingTrackCands[0];
               
	  for( unsigned j=1; j < overlappingTrackCands.size(); j++ ){
	    if( overlappingTrackCands[j]->getChi2Prob() > bestTrack->getChi2Prob() ){
	      delete bestTrack; //delete the old one, not needed anymore
	      bestTrack = overlappingTrackCands[j];
	    }
	    else{
	      delete overlappingTrackCands[j]; //delete this one
	    }
	  }
	  debug() << "Adding best track candidate with " << bestTrack->getHits().size() << " hits" << endmsg;
               
	  trackCandidates.push_back( bestTrack );
	}
      }
      else{ // we take all versions
	debug() << "Taking all " << overlappingTrackCands.size() << " versions of the track" << endmsg;
	trackCandidates.insert( trackCandidates.end(), overlappingTrackCands.begin(), overlappingTrackCands.end() );
      }
    }
      
    //if( _useCED ){
    //for( unsigned i=0; i < trackCandidates.size(); i++ ) KiTrackMarlin::drawTrackRandColor( trackCandidates[i] );
    //}
      
    /**********************************************************************************************/
    /*               Get the best subset of tracks                                                */
    /**********************************************************************************************/
      
    debug() << "The track candidates so far: " << endmsg;
    for( unsigned iTrack=0; iTrack < trackCandidates.size(); iTrack++ ){
      debug() << "track " << iTrack << ": " << trackCandidates[iTrack] << "\t" << KiTrackMarlin::getTrackHitInfo( trackCandidates[iTrack] ) << endmsg;
    }
      
    debug() << "\t\t---Get best subset of tracks---" << endmsg ;
      
    std::vector< ITrack* > tracks;
    std::vector< ITrack* > rejected;
      
    TrackCompatibilityShare1SP comp;
    TrackQIChi2Prob trackQI;
    TrackQISpecial trackQISpecial;
          
    if( _bestSubsetFinder == "SubsetHopfieldNN" ){
      debug() << "Use SubsetHopfieldNN for getting the best subset" << endmsg ;
         
      SubsetHopfieldNN< ITrack* > subset;
      subset.add( trackCandidates );
      subset.calculateBestSet( comp, trackQI );
      tracks = subset.getAccepted();
      rejected = subset.getRejected();
    }
    else if( _bestSubsetFinder == "SubsetSimple" ){
      debug() << "Use SubsetSimple for getting the best subset" << endmsg ;
         
      SubsetSimple< ITrack* > subset;
      subset.add( trackCandidates );
      subset.calculateBestSet( comp, trackQISpecial );
      tracks = subset.getAccepted();
      rejected = subset.getRejected();
    }
    else { // in any other case take all tracks
      debug() << "Input for subset = \"" << _bestSubsetFinder << "\". All tracks are kept" << endmsg ;
         
      tracks = trackCandidates;
    }
          
    if( _useCED ){
      //for( unsigned i=0; i < tracks.size(); i++ ) KiTrackMarlin::drawTrack( tracks[i] , 0x00ff00 );
      //for( unsigned i=0; i < rejected.size(); i++ ) KiTrackMarlin::drawTrack( rejected[i] , 0xff0000 );
    }
      
    for ( unsigned i=0; i<rejected.size(); i++){
      delete rejected[i];
    }
          
    /**********************************************************************************************/
    /*               Finally: Finalise and save the tracks                                        */
    /**********************************************************************************************/
    debug() << "\t\t---Save Tracks---" << endmsg ;
      
    //auto trkCol = _outColHdl.createAndPut();
    
    for (unsigned int i=0; i < tracks.size(); i++){
      FTDTrack* myTrack = dynamic_cast< FTDTrack* >( tracks[i] );
         
      if( myTrack != NULL ){
	edm4hep::Track trackImpl( *(myTrack->getLcioTrack()) );
            
	try{
	  finaliseTrack( &trackImpl );
	  //trkCol->addElement( trackImpl );
	  trkCol->push_back(trackImpl);
	}
	catch( FitterException& e ){
	  debug() << "ForwardTracking: track couldn't be finalized due to fitter error: " << e.what() << endmsg;
	  //delete trackImpl;
	}
      }
    }
    
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

    //evt->addCollection(trkCol,_ForwardTrackCollection.c_str());
          
    debug() << "Forward Tracking found and saved " << tracks.size() << " tracks in event " << _nEvt << endmsg; 
          
    /**********************************************************************************************/
    /*                Clean up                                                                    */
    /**********************************************************************************************/
    // delete all the created IHits
    for ( unsigned i=0; i<hitsTBD.size(); i++ )  delete hitsTBD[i];
      
    // delete the FTracks
    for (unsigned int i=0; i < tracks.size(); i++){ delete tracks[i];}
  }

  //if( _useCED ) MarlinCED::draw(this);
  
  _nEvt ++ ;
  return StatusCode::SUCCESS;
}

StatusCode ForwardTrackingAlg::finalize(){
  for ( unsigned i=0; i< _crit2Vec.size(); i++) delete _crit2Vec[i];
  for ( unsigned i=0; i< _crit3Vec.size(); i++) delete _crit3Vec[i];
  for ( unsigned i=0; i< _crit4Vec.size(); i++) delete _crit4Vec[i];
  _crit2Vec.clear();
  _crit3Vec.clear();
  _crit4Vec.clear();
   
  delete _sectorSystemFTD;
  _sectorSystemFTD = NULL;
  
  debug() << "There are " << _nTrackCandidates << "track candidates from CA and "<<  _nTrackCandidatesPlus
	  << " track Candidates with hits from overlapping hits" << endmsg
	  << "The ratio is " << float( _nTrackCandidatesPlus )/_nTrackCandidates << endmsg;

  return GaudiAlgorithm::finalize();
}

std::map< IHit* , std::vector< IHit* > > ForwardTrackingAlg::getOverlapConnectionMap(std::map< int , std::vector< IHit* > > & map_sector_hits, 
										  const SectorSystemFTD* secSysFTD,
										  float distMax){
      
  unsigned nConnections=0;
  
  std::map< IHit* , std::vector< IHit* > > map_hitFront_hitsBack;
  std::map< int , std::vector< IHit* > >::iterator it;
  
  //for every sector
  for ( it= map_sector_hits.begin() ; it != map_sector_hits.end(); it++ ){
    std::vector< IHit* > hitVecA = it->second;
    int sector = it->first;
    
    // get the neighbouring petals
    FTDNeighborPetalSecCon secCon( secSysFTD );
    std::set< int > targetSectors = secCon.getTargetSectors( sector );
          
    //for all neighbouring petals
    for ( std::set<int>::iterator itTarg = targetSectors.begin(); itTarg!=targetSectors.end(); itTarg++ ){
      std::vector< IHit* > hitVecB = map_sector_hits[ *itTarg ];
      for ( unsigned j=0; j < hitVecA.size(); j++ ){
	IHit* hitA = hitVecA[j];
	for ( unsigned k=0; k < hitVecB.size(); k++ ){
	  IHit* hitB = hitVecB[k];
                         
	  float dx = hitA->getX() - hitB->getX();
	  float dy = hitA->getY() - hitB->getY();
	  float dz = hitA->getZ() - hitB->getZ();
	  float dist = sqrt( dx*dx + dy*dy + dz*dz );
          
	  if (( dist < distMax )&& ( fabs( hitB->getZ() ) > fabs( hitA->getZ() ) )  ){ // if they are close enough and B is behind A
	    debug() << "Connected: (" << hitA->getX() << "," << hitA->getY() << "," << hitA->getZ() << ")-->("
		    << hitB->getX() << "," << hitB->getY() << "," << hitB->getZ() << ")" << endmsg;
	    
	    map_hitFront_hitsBack[ hitA ].push_back( hitB );
	    nConnections++;
	  }
	}
      } 
    }
  }
  debug() << "Connected " << map_hitFront_hitsBack.size() << " hits with " << nConnections << " possible overlapping hits" << endmsg;
      
  return map_hitFront_hitsBack;
}

std::string ForwardTrackingAlg::getInfo_map_sector_hits(){
  std::stringstream s;
  
  std::map< int , std::vector< IHit* > >::iterator it;
  
  for( it = _map_sector_hits.begin(); it != _map_sector_hits.end(); it++ ){
    std::vector<IHit*> hits = it->second;
    int sector = it->first;
    
    int side = _sectorSystemFTD->getSide( sector );
    unsigned layer = _sectorSystemFTD->getLayer( sector );
    unsigned module = _sectorSystemFTD->getModule( sector );
    unsigned sensor = _sectorSystemFTD->getSensor( sector );
    
    s << "sector " << sector  << " (si"
      << side << ",la"
      << layer << ",mo"
      << module << "se,"
      << sensor << ") has "
      << hits.size() << " hits\n" ;
  }  
     
  return s.str();   
}

std::vector < RawTrack > ForwardTrackingAlg::getRawTracksPlusOverlappingHits( RawTrack rawTrack , std::map< IHit* , std::vector< IHit* > >& map_hitFront_hitsBack ){
  // So we have a raw track (a vector of hits, that is) and a map, that tells us
  // for every hit, if there is another hit in the overlapping region behind it very close,
  // so that it could be part of the same track.
  //
  // We now want to find for a given track all possible tracks, when hits from the overlapping regions are added
  //
  // The method is this: start with pure track.
  // Make a vector of rawTracks and fill in the pure track.
  // For every hit on the original track do the following:
  // Check if there are overlapping hits.
  // For every overlapping hit take all the created tracks so far and make another version
  // with the overlapping hit added to it and add them to the vector of rawTracks.
  //
  //
  // Let's do an example: 
  // the original hits in the track are calles A,B and C.
  // A has one overlapping hit A1
  // and B has two overlapping hits B1 and B2.
  //
  // So we start with a vector containing only the original track: {(A,B,C)}
  //
  // We start with the first hit: A. It has one overlapping hit A1.
  // We take all tracks (which is the original one so far) and make another version containing A1 as well.
  // Then we add it to the vector of tracks:
  //
  // {(A,B,C)(A,A1,B,C)}
  //
  // On to the next hit from the original track: B. Here we have overlapping hits B1 and B2.
  // We take all the tracks so far and add versions with B1: (A,B,B1,C) and (A,A1,B,B1,C)
  // We don't immediately add them or otherwise, we would create a track containing B1 as well as B2, which is plainly wrong
  //
  // So instead we make the combinations with B2: (A,B,B2,C) and (A,A1,B,B2,C)
  // And now having gone through all overlapping hits of B, we add all the new versions to the vector:
  //
  // {(A,B,C)(A,A1,B,C)(A,B,B1,C)(A,A1,B,B1,C)(A,B,B2,C)(A,A1,B,B2,C)}
  //
  // So now we have all possible versions of the track with overlapping hits
   
  std::vector < RawTrack > rawTracksPlus;
   
  rawTracksPlus.push_back( rawTrack ); //add the original one
   
  // for every hit in the original track
  for( unsigned i=0; i < rawTrack.size(); i++ ){
    IHit* frontHit = rawTrack[i];
    
    // get the hits that are behind frontHit
    std::map< IHit* , std::vector< IHit* > >::iterator it;
    it = map_hitFront_hitsBack.find( frontHit );
    if( it == map_hitFront_hitsBack.end() ) continue; // if there are no hits on the back skip this one
    std::vector< IHit* > backHits = it->second; 
          
    // Create the different versions of the tracks so far with the hits from the back
      
    std::vector< RawTrack > newVersions; //here we store all the versions with overlapping hits from the different back hits at this frontHit
      
    // for every hit in back of the frontHit
    for( unsigned j=0; j<backHits.size(); j++ ){
      IHit* backHit = backHits[j];
      // for all tracks we have so far
      for( unsigned k=0; k<rawTracksPlus.size(); k++ ){
	RawTrack newVersion = rawTracksPlus[k];     // exact copy of the track
	newVersion.push_back( backHit );          // add the backHit to it   
	newVersions.push_back( newVersion );         // store it
        
      }
    }
    // Now put all the new versions of the tracks into the rawTracksPlus vector before we go on to the next
    // hit of the original track
    rawTracksPlus.insert( rawTracksPlus.end(), newVersions.begin(), newVersions.end() );
  }
     
  return rawTracksPlus;
}

bool ForwardTrackingAlg::setCriteria( unsigned round ){
  
  // delete the old ones
  for ( unsigned i=0; i< _crit2Vec.size(); i++) delete _crit2Vec[i];
  for ( unsigned i=0; i< _crit3Vec.size(); i++) delete _crit3Vec[i];
  for ( unsigned i=0; i< _crit4Vec.size(); i++) delete _crit4Vec[i];
  _crit2Vec.clear();
  _crit3Vec.clear();
  _crit4Vec.clear();
  
  bool newValuesGotUsed = false; // if new values are used
  for( unsigned i=0; i<_criteriaNames.size(); i++ ){
    std::string critName = _criteriaNames[i];

    float min = _critMinima[critName].back();
    float max = _critMaxima[critName].back();

    // use the value corresponding to the round, if there are no new ones for this criterion, just do nothing (the previous value stays in place)
    if( round + 1 <= _critMinima[critName].size() ){
      min =  _critMinima[critName][round];
      newValuesGotUsed = true;
    }
      
    if( round + 1 <= _critMaxima[critName].size() ){
      max =  _critMaxima[critName][round];
      newValuesGotUsed = true;
    }

    ICriterion* crit = Criteria::createCriterion( critName, min , max );

    // Some debug output about the created criterion
    std::string type = crit->getType();
    
    debug() <<  "Added: Criterion " << critName << " (type =  " << type 
	    << " ). Min = " << min
	    << ", Max = " << max
	    << ", round " << round << endmsg;
        
    // Add the new criterion to the corresponding vector
    if( type == "2Hit" ){
      _crit2Vec.push_back( crit );
    }
    else if( type == "3Hit" ){
      _crit3Vec.push_back( crit );
    }
    else if( type == "4Hit" ){
      _crit4Vec.push_back( crit );
    }
    else delete crit;
  }

  return newValuesGotUsed;
}

void ForwardTrackingAlg::finaliseTrack( edm4hep::Track* trackImpl ){
     
  Fitter fitter( trackImpl , _trkSystem );
   
  //trackImpl->trackStates().clear();
  int nState = trackImpl->trackStates_size();
  if(nState>0){
    debug() << "Track has " << nState << " TrackState now, should be cleared but not supported by EDM4hep" << endmsg;
    for(int i=0;i<nState;i++){
      debug() << trackImpl->getTrackStates(i).location << " ";
    }
    debug() << endmsg;
  }
  
  edm4hep::TrackState trkStateIP( *fitter.getTrackState( 1/*lcio::TrackState::AtIP*/ ) ) ;
  trkStateIP.location = 1;
  trackImpl->addToTrackStates( trkStateIP );
  
  edm4hep::TrackState trkStateFirstHit( *fitter.getTrackState( 2/*TrackState::AtFirstHit*/ ) ) ;
  trkStateFirstHit.location = 2;
  trackImpl->addToTrackStates( trkStateFirstHit );
  
  edm4hep::TrackState trkStateLastHit( *fitter.getTrackState( 3/*TrackState::AtLastHit*/ ) ) ;
  trkStateLastHit.location = 3;
  trackImpl->addToTrackStates( trkStateLastHit );
  
  edm4hep::TrackState trkStateAtCalo( *fitter.getTrackState( 4/*TrackState::AtCalorimeter*/ ) ) ;
  trkStateAtCalo.location = 4;
  trackImpl->addToTrackStates( trkStateAtCalo );
  
  trackImpl->setChi2( fitter.getChi2( 1 ) );
  trackImpl->setNdf(  fitter.getNdf ( 1 ) );
  
  const edm4hep::Vector3f& p = trkStateFirstHit.referencePoint;
  trackImpl->setRadiusOfInnermostHit( sqrt( p[0]*p[0] + p[1]*p[1] + p[2]*p[2] ) );
  
  std::map<int, int> hitNumbers; 
   
  hitNumbers[UTIL::ILDDetID::VXD] = 0;
  hitNumbers[UTIL::ILDDetID::SIT] = 0;
  hitNumbers[UTIL::ILDDetID::FTD] = 0;
  hitNumbers[UTIL::ILDDetID::TPC] = 0;
  hitNumbers[UTIL::ILDDetID::SET] = 0;
  hitNumbers[UTIL::ILDDetID::ETD] = 0;
  
  unsigned int nHits = trackImpl->trackerHits_size();
  for( unsigned j=0; j<nHits; j++ ){
    const edm4hep::ConstTrackerHit& hit = trackImpl->getTrackerHits(j);
    UTIL::BitField64 encoder( UTIL::ILDCellID0::encoder_string );
    encoder.setValue( hit.getCellID() );
    int subdet =  encoder[UTIL::ILDCellID0::subdet];
        
    ++hitNumbers[ subdet ];
  }
  
  //trackImpl->subdetectorHitNumbers().resize(2 * lcio::ILDDetID::ETD);
  //trackImpl->subdetectorHitNumbers()[ 2 * lcio::ILDDetID::VXD - 2 ] = hitNumbers[lcio::ILDDetID::VXD];
  //trackImpl->subdetectorHitNumbers()[ 2 * lcio::ILDDetID::FTD - 2 ] = hitNumbers[lcio::ILDDetID::FTD];
  //trackImpl->subdetectorHitNumbers()[ 2 * lcio::ILDDetID::SIT - 2 ] = hitNumbers[lcio::ILDDetID::SIT];
  //trackImpl->subdetectorHitNumbers()[ 2 * lcio::ILDDetID::TPC - 2 ] = hitNumbers[lcio::ILDDetID::TPC];
  //trackImpl->subdetectorHitNumbers()[ 2 * lcio::ILDDetID::SET - 2 ] = hitNumbers[lcio::ILDDetID::SET];
  //trackImpl->subdetectorHitNumbers()[ 2 * lcio::ILDDetID::ETD - 2 ] = hitNumbers[lcio::ILDDetID::ETD];
  //trackImpl->subdetectorHitNumbers()[ 2 * lcio::ILDDetID::VXD - 1 ] = hitNumbers[lcio::ILDDetID::VXD];
  //trackImpl->subdetectorHitNumbers()[ 2 * lcio::ILDDetID::FTD - 1 ] = hitNumbers[lcio::ILDDetID::FTD];
  //trackImpl->subdetectorHitNumbers()[ 2 * lcio::ILDDetID::SIT - 1 ] = hitNumbers[lcio::ILDDetID::SIT];
  //trackImpl->subdetectorHitNumbers()[ 2 * lcio::ILDDetID::TPC - 1 ] = hitNumbers[lcio::ILDDetID::TPC];
  //trackImpl->subdetectorHitNumbers()[ 2 * lcio::ILDDetID::SET - 1 ] = hitNumbers[lcio::ILDDetID::SET];
  //trackImpl->subdetectorHitNumbers()[ 2 * lcio::ILDDetID::ETD - 1 ] = hitNumbers[lcio::ILDDetID::ETD];
  trackImpl->addToSubDetectorHitNumbers(hitNumbers[UTIL::ILDDetID::VXD]);
  trackImpl->addToSubDetectorHitNumbers(hitNumbers[UTIL::ILDDetID::SIT]);
  trackImpl->addToSubDetectorHitNumbers(hitNumbers[UTIL::ILDDetID::FTD]);
  trackImpl->addToSubDetectorHitNumbers(hitNumbers[UTIL::ILDDetID::TPC]);
  trackImpl->addToSubDetectorHitNumbers(hitNumbers[UTIL::ILDDetID::SET]);
  trackImpl->addToSubDetectorHitNumbers(hitNumbers[UTIL::ILDDetID::ETD]);
     
  return;
}


