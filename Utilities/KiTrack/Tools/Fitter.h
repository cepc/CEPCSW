#ifndef Fitter_h
#define Fitter_h

#include "TrackSystemSvc/IMarlinTrkSystem.h"
#include "TrackSystemSvc/IMarlinTrack.h"
#include "edm4hep/Track.h"
//#include "lcio.h"

#include "Math/ProbFunc.h"


//using namespace lcio;



class FitterException : public std::exception {
   
   
protected:
   std::string message ;
   
   FitterException(){  /*no_op*/ ; } 
   
public: 
   virtual ~FitterException() { /*no_op*/; } 
   
   FitterException( const std::string& text ){
      message = "FitterException: " + text ;
   }
   
   virtual const char* what() const noexcept { return  message.c_str() ; } 
   
};



/** A class to store additional information together with a TrackState, namely the chi2 value and Ndf
 */
class TrackStatePlus{
   
   
   
public: 
   
  const edm4hep::TrackState* getTrackState() const {return _trackState; }
  double getChi2() const {return _chi2;}
  int getNdf() const {return _Ndf;}
  
 TrackStatePlus( const edm4hep::TrackState* trackState, double chi2, int Ndf ):
  _trackState( trackState ), _chi2( chi2 ), _Ndf( Ndf ){}
   
   
   
private: 
   
  const edm4hep::TrackState* _trackState;
   double _chi2;
   int _Ndf;
   
   
};


/** A class to make it quick to fit a track or hits and get back the chi2, Ndf and chi2prob values and
 * also bundle the code used for that, so it doesn't have to be copied all over the places.
 */
class Fitter{
   
   
public:
   
  Fitter( edm4hep::Track* track , MarlinTrk::IMarlinTrkSystem* trkSystem );
  Fitter( std::vector < edm4hep::ConstTrackerHit > trackerHits, MarlinTrk::IMarlinTrkSystem* trkSystem );
  Fitter( edm4hep::Track* track , MarlinTrk::IMarlinTrkSystem* trkSystem, int VXDFlag );  

   
   double getChi2Prob( int trackStateLocation ) ;
   double getChi2( int trackStateLocation ) ;
   int getNdf( int trackStateLocation ) ;
   
   //TODO: maybe add methods for custom points (location: TrackState::AtOther) In that case, the point would have to
   // be passed as well. (or only the point is passed)
   
   const edm4hep::TrackState* getTrackState( int trackStateLocation ) ;
   
   ~Fitter(){ 
      
      for( unsigned i=0; i<_trackStatesPlus.size(); i++ ){
         
         delete _trackStatesPlus[i]->getTrackState();
         delete _trackStatesPlus[i];
         
      }
      _trackStatesPlus.clear();
      
      delete _marlinTrk;
      
   }
   
private:
   
   void init_BField();

   const TrackStatePlus* getTrackStatePlus( int trackStateLocation ) ;

   void fit();

   void fitVXD();
   
   static float _bField;
   
   
   std::vector< edm4hep::ConstTrackerHit > _trackerHits;
   
   /** here the created TrackStates (plus) are stored */
   std::vector< const TrackStatePlus* > _trackStatesPlus;
   
   MarlinTrk::IMarlinTrkSystem* _trkSystem;
   
   MarlinTrk::IMarlinTrack* _marlinTrk;
   
   // No copy constructor or assignment needed so far. so private for safety
   // If they are needed, they need to be implemented in a clean way first!
   Fitter( const Fitter& f ){};
   Fitter& operator= ( Fitter const& f ){return *this;}
   
};

#endif
