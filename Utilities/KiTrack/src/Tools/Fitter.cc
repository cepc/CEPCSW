#include "Tools/Fitter.h"

#include <algorithm>

//#include "marlin/Global.h"
//#include "UTIL/LCTrackerConf.h"
#include <UTIL/ILDConf.h>

#include <DD4hep/DD4hepUnits.h>
#include <DD4hep/Detector.h>

#include "TrackSystemSvc/HelixTrack.h"
#include "DataHelper/Navigation.h"
#include "Tools/KiTrackMarlinTools.h"

typedef std::vector<edm4hep::ConstTrackerHit> TrackerHitVec;
using namespace MarlinTrk;

// by fucd: 3.5->3.0 default, will be read from GeomSvc
// if compare to Marlin, should change to 3.5
float Fitter::_bField = 3.5;//later on overwritten with the value read by geo file

void Fitter::init_BField(){

  // B field from DD4hep
  /*
  dd4hep::Detector & lcdd = dd4hep::Detector::getInstance();
  const double pos[3]={0,0,0}; 
  double bFieldVec[3]={0,0,0}; 
  lcdd.field().magneticField(pos,bFieldVec); // get the magnetic field vector from DD4hep
  _bField = bFieldVec[2]/dd4hep::tesla; // z component at (0,0,0)    
  */

}

bool compare_TrackerHit_z( edm4hep::ConstTrackerHit a, edm4hep::ConstTrackerHit b ){
  return ( fabs(a.getPosition()[2]) < fabs( b.getPosition()[2]) ); //compare their z values
}

bool compare_TrackerHit_R( edm4hep::ConstTrackerHit a, edm4hep::ConstTrackerHit b ){
  double Rad_a2 = (a.getPosition()[0]*a.getPosition()[0]) + (a.getPosition()[1]*a.getPosition()[1]) ;
  double Rad_b2 = (b.getPosition()[0]*b.getPosition()[0]) + (b.getPosition()[1]*b.getPosition()[1]) ;
  
  return ( Rad_a2 < Rad_b2 ); //compare their radii
}

Fitter::Fitter( edm4hep::Track* track , MarlinTrk::IMarlinTrkSystem* trkSystem ): _trkSystem( trkSystem ){
  _trackerHits.clear();

  std::copy(track->trackerHits_begin(), track->trackerHits_end(), std::back_inserter(_trackerHits));
  //_trackerHits = track->getTrackerHits();

  fit();
}

Fitter::Fitter( edm4hep::Track* track , MarlinTrk::IMarlinTrkSystem* trkSystem, int VXDFlag ): _trkSystem( trkSystem ){
  _trackerHits.clear();
  std::copy(track->trackerHits_begin(), track->trackerHits_end(), std::back_inserter(_trackerHits));
  //_trackerHits = track->getTrackerHits();
  fitVXD();
}

Fitter::Fitter( std::vector<edm4hep::ConstTrackerHit> trackerHits , MarlinTrk::IMarlinTrkSystem* trkSystem ): _trkSystem( trkSystem ){
  _trackerHits = trackerHits;
  fit();
}

void Fitter::fitVXD(){
  //create the MarlinTrk
  _marlinTrk = _trkSystem->createTrack();

  /**********************************************************************************************/
  /*       Add the hits to the MarlinTrack                                                      */
  /**********************************************************************************************/
   
  // hits are in reverse order 
  std::sort( _trackerHits.begin(), _trackerHits.end(), KiTrackMarlin::compare_TrackerHit_R );
  // now at [0] is the hit with the smallest |z| and at [1] is the one with a bigger |z| and so on
  // So the direction of the hits when following the index from 0 on is:
  // from inside out: from the IP into the distance. 
  // (It is important to keep in mind, in which direction we fit, when using MarlinTrk)
  
  TrackerHitVec::iterator it;
     
  unsigned number_of_added_hits = 0;
  unsigned ndof_added = 0;
  std::vector< edm4hep::ConstTrackerHit > added_hits;
  std::vector< edm4hep::ConstTrackerHit > added_hits_2D;
  
  for( it = _trackerHits.begin() ; it != _trackerHits.end() ; ++it ) {
    edm4hep::ConstTrackerHit trkHit = Navigation::Instance()->GetTrackerHit((*it).getObjectID());
    bool isSuccessful = false; 
    
    if( UTIL::BitSet32( trkHit.getType() )[ UTIL::ILDTrkHitTypeBit::COMPOSITE_SPACEPOINT ]   ){ //it is a composite spacepoint
      //Split it up and hits to the MarlinTrk
      std::vector< edm4hep::ConstTrackerHit > rawHits;
      //const LCObjectVec rawObjects = trkHit.getRawHits();
      //for( unsigned k=0; k<rawObjects.size(); k++ ) rawHits.push_back( dynamic_cast< ConstTrackerHit >( rawObjects[k] ) );
      int nRawHit = trkHit.rawHits_size();
      for( unsigned k=0; k< nRawHit; k++ ){
	edm4hep::ConstTrackerHit rawHit = Navigation::Instance()->GetTrackerHit(trkHit.getRawHits(k));
	rawHits.push_back(rawHit);
      }
      std::sort( rawHits.begin(), rawHits.end(), compare_TrackerHit_R );
      
      for( unsigned k=0; k< rawHits.size(); k++ ){
	if( _marlinTrk->addHit( rawHits[k] ) == IMarlinTrack::success ){
	  isSuccessful = true; //if at least one hit from the spacepoint gets added
	  ++ndof_added; // 1 degree of freedom for each strip hit
	}
	else{
	  //std::cout << "Cannot addHit " << rawHits[k]->id() << " to MarlinTrk" << std::endl; 
	}
      }
    }
    else { // normal non composite hit
      if (_marlinTrk->addHit( trkHit ) == 0) {
	isSuccessful = true;
	ndof_added += 2;
      }
    }
    
    if (isSuccessful) {
      added_hits.push_back(trkHit);
      ++number_of_added_hits;
    }
    else{
      //std::cout << "DEBUG Fitter::fit(): Hit " << it - _trackerHits.begin() << " Dropped " << std::endl;          
    }
  }
  
  if( ndof_added < 6 ) {
    std::stringstream s;
    s << "Fitter::fit(): Cannot fit less with less than 6 degrees of freedom. Number of hits =  " << number_of_added_hits << " ndof = " << ndof_added << "\n";
    
    throw FitterException( s.str() );
  }
   
  /**********************************************************************************************/
  /*       Create a helix from the first, last and middle hit                                   */
  /**********************************************************************************************/
  
  for (unsigned ihit=0; ihit <added_hits.size(); ++ihit) {
    // check if this a space point or 2D hit 
    if(UTIL::BitSet32( added_hits[ihit].getType() )[ UTIL::ILDTrkHitTypeBit::ONE_DIMENSIONAL ] == false ){
      // then add to the list 
      added_hits_2D.push_back(added_hits[ihit]);
      
    }
  }
  
  // initialise with space-points not strips 
  // make a helix from 3 hits to get a trackstate
  const edm4hep::Vector3d x1 = added_hits_2D[0].getPosition();
  const edm4hep::Vector3d x2 = added_hits_2D[ added_hits_2D.size()/2 ].getPosition();
  const edm4hep::Vector3d x3 = added_hits_2D.back().getPosition();   
  
  init_BField();
  HelixTrack helixTrack( x1, x2, x3, _bField, HelixTrack::forwards );
  
  helixTrack.moveRefPoint(0.0, 0.0, 0.0);
  
  //const float referencePoint[3] = { float(helixTrack.getRefPointX()) , float(helixTrack.getRefPointY() ), float(helixTrack.getRefPointZ() )};
  edm4hep::Vector3f referencePoint = { float(helixTrack.getRefPointX()) , float(helixTrack.getRefPointY() ), float(helixTrack.getRefPointZ() )};
  
  /**********************************************************************************************/
  /*       Create a TrackStateImpl from the helix values and use it to initalise the fit        */
  /**********************************************************************************************/
  std::array<float,15> covMatrix;
  
  for (unsigned icov = 0; icov<covMatrix.size(); ++icov) {
    covMatrix[icov] = 0;
  }
  
  covMatrix[0]  = ( 1.e6 ); //sigma_d0^2
  covMatrix[2]  = ( 1.e2 ); //sigma_phi0^2
  covMatrix[5]  = ( 1.e-4 ); //sigma_omega^2
  covMatrix[9]  = ( 1.e6 ); //sigma_z0^2
  covMatrix[14] = ( 1.e2 ); //sigma_tanl^2
  
  
  edm4hep::TrackState trackState = {0/*TrackState::AtOther*/, 
				    helixTrack.getD0(), 
				    helixTrack.getPhi0(), 
				    helixTrack.getOmega(), 
				    helixTrack.getZ0(), 
				    helixTrack.getTanLambda(),
				    referencePoint,
				    covMatrix};
  
  //init_BField();
  _marlinTrk->initialise( trackState, _bField, IMarlinTrack::backward ) ;
  
  //     _marlinTrk->initialise( IMarlinTrack::backward ) ;
   
  /**********************************************************************************************/
  /*       Do the fit                                                                           */
  /**********************************************************************************************/
  
  int fit_status = 0;
  
  try{
    
    fit_status = _marlinTrk->fit() ; 
    
  }
  catch( MarlinTrk::Exception& e ){
    
    std::stringstream s;
    s << "Fitter::fit(): Couldn't fit, MarlinTrk->fit() gave: " << e.what() << "\n";
    throw FitterException( s.str() );
    
  }
  
  if( fit_status != IMarlinTrack::success ){ 
    
    std::stringstream s;
    s << "Fitter::fit(): MarlinTrk->fit() wasn't successful, fit_status = " << fit_status << "\n";
    throw FitterException( s.str() );
    
  }
  
  
  // fitting finished get hits in the fit for safety checks:
  
  std::vector<std::pair<edm4hep::ConstTrackerHit, double> > hits_in_fit;
  
  // remember the hits are ordered in the order in which they were fitted
  // here we are fitting inwards so the first is the last and vice verse
  
  _marlinTrk->getHitsInFit(hits_in_fit);
  
  if( hits_in_fit.size() < 3 ) {
    
    
    std::stringstream s;
    s << "Fitter::fit() Less than 3 hits in fit: Only " << hits_in_fit.size() << 
      " of " << _trackerHits.size() << " hits\n";
    
    throw FitterException( s.str() );
    
  }
  edm4hep::ConstTrackerHit first_hit_in_fit = hits_in_fit.back().first;
  if (! first_hit_in_fit.isAvailable()) {
    throw FitterException( std::string("Fitter::fit(): TrackerHit pointer to first hit == NULL ")  ) ;
  }
  
  
  edm4hep::ConstTrackerHit last_hit_in_fit = hits_in_fit.front().first;
  if (!last_hit_in_fit.isAvailable()) {
    throw FitterException( std::string("Fitter::fit(): TrackerHit pointer to last hit == NULL ")  ) ;
  }
  
  return;
}



void Fitter::fit(){
  //create the MarlinTrk
  _marlinTrk = _trkSystem->createTrack();
  
  /**********************************************************************************************/
  /*       Add the hits to the MarlinTrack                                                      */
  /**********************************************************************************************/
  
  // hits are in reverse order 
  std::sort( _trackerHits.begin(), _trackerHits.end(), KiTrackMarlin::compare_TrackerHit_z );
  // now at [0] is the hit with the smallest |z| and at [1] is the one with a bigger |z| and so on
  // So the direction of the hits when following the index from 0 on is:
  // from inside out: from the IP into the distance. 
  // (It is important to keep in mind, in which direction we fit, when using MarlinTrk)
  
  TrackerHitVec::iterator it;
  
  unsigned number_of_added_hits = 0;
  unsigned ndof_added = 0;
  std::vector<edm4hep::ConstTrackerHit> added_hits;
  
  for( it = _trackerHits.begin() ; it != _trackerHits.end() ; ++it ) {
    edm4hep::ConstTrackerHit trkHit = Navigation::Instance()->GetTrackerHit((*it).getObjectID());
    bool isSuccessful = false; 
    //std::cout << "Hit " << trkHit->id() << " " << trkHit.getPosition() << std::endl;
    if( UTIL::BitSet32( trkHit.getType() )[ UTIL::ILDTrkHitTypeBit::COMPOSITE_SPACEPOINT ]   ){ //it is a composite spacepoint
      //Split it up and hits to the MarlinTrk
      std::vector<edm4hep::ConstTrackerHit> rawHits;
      //const LCObjectVec rawObjects = trkHit.getRawHits();                    
      //for( unsigned k=0; k<rawObjects.size(); k++ ) rawHits.push_back( dynamic_cast< ConstTrackerHit >( rawObjects[k] ) );
      int nRawHit = trkHit.rawHits_size();
      for( unsigned k=0; k< nRawHit; k++ ){
	edm4hep::ConstTrackerHit rawHit = Navigation::Instance()->GetTrackerHit(trkHit.getRawHits(k));
	//std::cout << "Raw Hit " << rawHit->id() << " " << rawHit.getPosition() << std::endl;
	rawHits.push_back(rawHit);
      }
      std::sort( rawHits.begin(), rawHits.end(), compare_TrackerHit_z );
      
      for( unsigned k=0; k< rawHits.size(); k++ ){
	if( _marlinTrk->addHit( rawHits[k] ) == IMarlinTrack::success ){
	   isSuccessful = true; //if at least one hit from the spacepoint gets added
	   ++ndof_added; // 1 degree of freedom for each strip hit
	}
	else{
	  //std::cout << "Cannot addHit " << rawHits[k]->id() << " to MarlinTrk" << std::endl;
	}
      }
    }
    else { // normal non composite hit
      
      if (_marlinTrk->addHit( trkHit ) == 0) {
	isSuccessful = true;
	ndof_added += 2;
      }
    }
    
    if (isSuccessful) {
      added_hits.push_back(trkHit);
      ++number_of_added_hits;
    }
    else{
      //std::cout << "DEBUG: Fitter::fit(): Hit " << it - _trackerHits.begin() << " Dropped " << std::endl;          
    }
  }
  
  if( ndof_added < 6 ) {
    std::stringstream s;
    s << "Fitter::fit(): Cannot fit less with less than 6 degrees of freedom. Number of hits =  " << number_of_added_hits << " ndof = " << ndof_added << "\n";
    
    throw FitterException( s.str() );
  }
   
  /**********************************************************************************************/
  /*       Create a helix from the first, last and middle hit                                   */
  /**********************************************************************************************/
    
  // initialise with space-points not strips 
  // make a helix from 3 hits to get a trackstate
  const edm4hep::Vector3d x1 = added_hits[0].getPosition();
  const edm4hep::Vector3d x2 = added_hits[ added_hits.size()/2 ].getPosition();
  const edm4hep::Vector3d x3 = added_hits.back().getPosition();
  
  init_BField();
  HelixTrack helixTrack( x1, x2, x3, _bField, HelixTrack::forwards );
  
  helixTrack.moveRefPoint(0.0, 0.0, 0.0);
  
  //const float referencePoint[3] = { float(helixTrack.getRefPointX()) , float(helixTrack.getRefPointY()) , float(helixTrack.getRefPointZ()) };
  edm4hep::Vector3f referencePoint = { float(helixTrack.getRefPointX()) , float(helixTrack.getRefPointY()) , float(helixTrack.getRefPointZ()) };
  
  /**********************************************************************************************/
  /*       Create a TrackStateImpl from the helix values and use it to initalise the fit        */
  /**********************************************************************************************/
  std::array<float,15> covMatrix;
  for (unsigned icov = 0; icov<covMatrix.size(); ++icov) {
    covMatrix[icov] = 0;
  }
  
  covMatrix[0]  = ( 1.e6 ); //sigma_d0^2
  covMatrix[2]  = ( 1.e2 ); //sigma_phi0^2
  covMatrix[5]  = ( 1.e-4 ); //sigma_omega^2
  covMatrix[9]  = ( 1.e6 ); //sigma_z0^2
  covMatrix[14] = ( 1.e2 ); //sigma_tanl^2
  
  edm4hep::TrackState trackState = {0/*TrackState::AtOther*/, 
				    helixTrack.getD0(), 
				    helixTrack.getPhi0(), 
				    helixTrack.getOmega(), 
				    helixTrack.getZ0(), 
				    helixTrack.getTanLambda(), 
				    referencePoint,
				    covMatrix};
  
  //init_BField();
  _marlinTrk->initialise( trackState, _bField, IMarlinTrack::backward ) ;
  
  //     _marlinTrk->initialise( IMarlinTrack::backward ) ;
  
  /**********************************************************************************************/
  /*       Do the fit                                                                           */
  /**********************************************************************************************/
  int fit_status = 0;
  
  try{
    fit_status = _marlinTrk->fit() ; 
  }
  catch( MarlinTrk::Exception& e ){
    std::stringstream s;
    s << "Fitter::fit(): Couldn't fit, MarlinTrk->fit() gave: " << e.what() << "\n";
    throw FitterException( s.str() );
  }
   
  if( fit_status != IMarlinTrack::success ){ 
    std::stringstream s;
    s << "Fitter::fit(): MarlinTrk->fit() wasn't successful, fit_status = " << fit_status << "\n";
    throw FitterException( s.str() );
  }
  
  // fitting finished get hits in the fit for safety checks:
  
  std::vector<std::pair<edm4hep::ConstTrackerHit, double> > hits_in_fit;
  
  // remember the hits are ordered in the order in which they were fitted
  // here we are fitting inwards so the first is the last and vice verse
  
  _marlinTrk->getHitsInFit(hits_in_fit);
  
  if( hits_in_fit.size() < 3 ) {
    std::stringstream s;
    s << "Fitter::fit() Less than 3 hits in fit: Only " << hits_in_fit.size()
      << " of " << _trackerHits.size() << " hits\n";
    
    throw FitterException( s.str() );
  }
  edm4hep::ConstTrackerHit first_hit_in_fit = hits_in_fit.back().first;
  if (!first_hit_in_fit.isAvailable()) {
    throw FitterException( std::string("Fitter::fit(): TrackerHit pointer to first hit == NULL ")  ) ;
  }
  
  edm4hep::ConstTrackerHit last_hit_in_fit = hits_in_fit.front().first;
  if (!last_hit_in_fit.isAvailable()) {
    throw FitterException( std::string("Fitter::fit(): TrackerHit pointer to last hit == NULL ")  ) ;
  }
  
  return;
}
   
   
const edm4hep::TrackState* Fitter::getTrackState( int trackStateLocation ){
  return getTrackStatePlus( trackStateLocation )->getTrackState();
}

double Fitter::getChi2Prob( int trackStateLocation ){
  return ROOT::Math::chisquared_cdf_c( getChi2( trackStateLocation ) , getNdf( trackStateLocation ) );   
}

double Fitter::getChi2( int trackStateLocation ){
  return getTrackStatePlus( trackStateLocation )->getChi2();
}
   
int Fitter::getNdf( int trackStateLocation ){
   return getTrackStatePlus( trackStateLocation )->getNdf();
}
   
   
const TrackStatePlus* Fitter::getTrackStatePlus( int trackStateLocation ){
   // check if there is already an entry with this trackState location
   for( unsigned i=0; i<_trackStatesPlus.size(); i++ ){
     if( _trackStatesPlus[i]->getTrackState()->location == trackStateLocation ){
       return _trackStatesPlus[i];
     }
   }
   
   // If we reach this point, obviously no trackState with the given location has been created so far
   // Thus we create it now
   edm4hep::TrackState* trackState = new edm4hep::TrackState;
   int return_code = 0;
   double chi2;
   int ndf;
   switch( trackStateLocation ){
   case 1/*lcio::TrackState::AtIP*/:{
     const edm4hep::Vector3d point(0.,0.,0.); // nominal IP
              
     return_code = _marlinTrk->propagate(point, *trackState, chi2, ndf ) ;
         
     if (return_code != MarlinTrk::IMarlinTrack::success ) {
       delete trackState;
            
       std::stringstream s;
       s << "Fitter::getTrackStatePlus(): Couldn't create TrackState at IP, return code from propagation = " << return_code << "\n";
       throw FitterException( s.str() );
       
       break;
     }
     else{
       trackState->location = trackStateLocation;
       TrackStatePlus* trackStatePlus = new TrackStatePlus( trackState, chi2, ndf );
       _trackStatesPlus.push_back( trackStatePlus );
       return trackStatePlus;
     }
   }
   case 2/*lcio::TrackState::AtFirstHit*/:{
     std::vector<std::pair<edm4hep::ConstTrackerHit, double> > hits_in_fit;
         
     // remember the hits are ordered in the order in which they were fitted
     // here we are fitting inwards so the first is the last and vice verse
     _marlinTrk->getHitsInFit(hits_in_fit);
     
     edm4hep::ConstTrackerHit first_hit_in_fit = hits_in_fit.back().first;
          
     return_code = _marlinTrk->getTrackState(first_hit_in_fit, *trackState, chi2, ndf ) ;
     
     if(return_code !=MarlinTrk::IMarlinTrack::success){
       
       delete trackState;
       
       std::stringstream s;
       s << "Fitter::getTrackStatePlus(): Couldn't create TrackState at first hit, return code from propagation = " << return_code << "\n";
       throw FitterException( s.str() );
	    
       break;
     }
     else{
       trackState->location = trackStateLocation;
       TrackStatePlus* trackStatePlus = new TrackStatePlus( trackState, chi2, ndf );
       _trackStatesPlus.push_back( trackStatePlus );
       return trackStatePlus;
     }
   }
   case 3/*lcio::TrackState::AtLastHit*/:{
     std::vector<std::pair<edm4hep::ConstTrackerHit, double> > hits_in_fit;
     _marlinTrk->getHitsInFit(hits_in_fit);
     
     edm4hep::ConstTrackerHit last_hit_in_fit = hits_in_fit.front().first;
          
     return_code = _marlinTrk->getTrackState(last_hit_in_fit, *trackState, chi2, ndf ) ;
         
     if(return_code !=MarlinTrk::IMarlinTrack::success){
       delete trackState;
            
       std::stringstream s;
       s << "Fitter::getTrackStatePlus(): Couldn't create TrackState at last hit, return code from propagation = " << return_code << "\n";
       throw FitterException( s.str() );
       
       break;
     }
     else{
       trackState->location = trackStateLocation;
       TrackStatePlus* trackStatePlus = new TrackStatePlus( trackState, chi2, ndf );
       _trackStatesPlus.push_back( trackStatePlus );
       return trackStatePlus;
     }
     break;
   }
   case 4/*lcio::TrackState::AtCalorimeter*/:{
     std::vector<std::pair<edm4hep::ConstTrackerHit, double> > hits_in_fit;
     _marlinTrk->getHitsInFit(hits_in_fit);
     
     edm4hep::ConstTrackerHit last_hit_in_fit = hits_in_fit.front().first;
          
     UTIL::BitField64 encoder( UTIL::ILDCellID0::encoder_string ) ; 
     encoder.reset() ;  // reset to 0
         
     // ================== need to get the correct ID(s) for the calorimeter face  ============================
     unsigned ecal_barrel_face_ID = UTIL::ILDDetID::ECAL ;
     //unsigned ecal_endcap_face_ID = UTIL::ILDDetID::ECAL_ENDCAP ;
     //=========================================================================================================

     encoder[UTIL::ILDCellID0::subdet] =  ecal_barrel_face_ID ;
     encoder[UTIL::ILDCellID0::side]   = UTIL::ILDDetID::barrel;
     encoder[UTIL::ILDCellID0::layer]  = 0 ;
     
     int detElementID = 0;
     return_code = _marlinTrk->propagateToLayer(encoder.lowWord(), last_hit_in_fit, *trackState, chi2, ndf, detElementID, IMarlinTrack::modeForward ) ;
     
     if (return_code == MarlinTrk::IMarlinTrack::no_intersection ) { // try forward or backward
       //encoder[UTIL::ILDCellID0::subdet] = ecal_endcap_face_ID ;
       
       const edm4hep::TrackState* trkStateLastHit = getTrackStatePlus( 3/*lcio::TrackState::AtLastHit*/ )->getTrackState();
       
       if (trkStateLastHit->tanLambda>0) {
	 encoder[UTIL::ILDCellID0::side] = UTIL::ILDDetID::fwd;
       }
       else{
	 encoder[UTIL::ILDCellID0::side] = UTIL::ILDDetID::bwd;
       }
       return_code = _marlinTrk->propagateToLayer(encoder.lowWord(), last_hit_in_fit, *trackState, chi2, ndf, detElementID, IMarlinTrack::modeForward ) ;
     }
                 
     if(return_code !=MarlinTrk::IMarlinTrack::success){
       delete trackState;
            
       std::stringstream s;
       s << "Fitter::getTrackStatePlus(): Couldn't create TrackState at Calorimeter, return code from propagation = " << return_code << "\n";
       throw FitterException( s.str() );
       break;
     }
     else{
       trackState->location = trackStateLocation;
       TrackStatePlus* trackStatePlus = new TrackStatePlus( trackState, chi2, ndf );
       _trackStatesPlus.push_back( trackStatePlus );
       return trackStatePlus;
     }
   }
   default:{
     std::stringstream s;
     s << "Creation of a trackState for the given location " << trackStateLocation 
       << " is not yet implemented for the class Fitter. \nImplemented are: AtIP, AtFirstHit, AtLastHit, AtCalorimeter.\n"
       << "If another location is desired, it must be implemented in the method Fitter::getTrackStatePlus.\n";
     
     throw FitterException( s.str() );
     return NULL;
   }
   }
      
   return NULL; // if we haven't returned so far, there was no success, so we return NULL
}
   
   

