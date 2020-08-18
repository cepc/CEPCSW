#include "ILDImpl/FTDTrack.h"


#include <algorithm>
#include <math.h>

//#include "UTIL/LCTrackerConf.h"

// Root, for calculating the chi2 probability. 
#include "Math/ProbFunc.h"



using namespace KiTrackMarlin;

/** @return if the absolute z value of hit a is bigger than that of hit b */
bool compare_IHit_z( IHit* a, IHit* b ){
   
   return ( fabs( a->getZ() ) < fabs( b->getZ() ) ); //compare their z values
   
}



FTDTrack::FTDTrack( MarlinTrk::IMarlinTrkSystem* trkSystem ){
   
   _trkSystem = trkSystem;
   _chi2Prob = 0.;
 
   _lcioTrack = new edm4hep::Track();
   
   
}

FTDTrack::FTDTrack( std::vector< IFTDHit* > hits , MarlinTrk::IMarlinTrkSystem* trkSystem ){
   
   
   _trkSystem = trkSystem;
   _chi2Prob = 0.;
   
   _lcioTrack = new edm4hep::Track();
   
   for( unsigned i=0; i < hits.size(); i++ ){
      
      addHit( hits[i] );
      
      
   }
   
}


FTDTrack::FTDTrack( const FTDTrack& f ){

   //make a new copied lcio track
  _lcioTrack = new edm4hep::Track( *f._lcioTrack );
   
   
   _hits = f._hits;
   _chi2Prob = f._chi2Prob;
   _trkSystem = f._trkSystem;

}

FTDTrack & FTDTrack::operator= (const FTDTrack & f){
   
   if (this == &f) return *this;   //protect against self assignment
   
   //make a new copied lcio track
   _lcioTrack = new edm4hep::Track( *f._lcioTrack );
   
   
   _hits = f._hits;
   _chi2Prob = f._chi2Prob;
   _trkSystem = f._trkSystem;
   
   return *this;
   
}



void FTDTrack::addHit( IFTDHit* hit ){
  if ( hit != NULL ){
    _hits.push_back( hit );
    // and sort the track again
    sort( _hits.begin(), _hits.end(), compare_IHit_z );
    _lcioTrack->addToTrackerHits( *hit->getTrackerHit() );
  }
}





void FTDTrack::fit() {
   
   
   Fitter fitter( _lcioTrack , _trkSystem );
   
   
   _lcioTrack->setChi2( fitter.getChi2( 1/*by fucd AtIP=1 in LCIO, changed to CepC rule in future: lcio::TrackState::AtIP*/ ) );
   _lcioTrack->setNdf( fitter.getNdf( 1/*lcio::TrackState::AtIP*/ ) );
   _chi2Prob = fitter.getChi2Prob( 1/*lcio::TrackState::AtIP*/ );
   
   edm4hep::TrackState trkState( *fitter.getTrackState( 1/*lcio::TrackState::AtIP*/ ) ) ;
   trkState.location = 1/*lcio::TrackState::AtIP*/ ;
   _lcioTrack->addToTrackStates( trkState );
   
   
}


double FTDTrack::getQI() const{
  
   
   double QI = _chi2Prob;
   
   // make sure QI is between 0 and 1
   if (QI > 1. ) QI = 1.;
   if (QI < 0. ) QI = 0.;
   
   return QI;
   
}







