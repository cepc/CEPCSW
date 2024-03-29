#ifndef VXDHelixFitter_h
#define VXDHelixFitter_h

#include "edm4hep/Track.h"
#include "edm4hep/TrackerHit.h"

//#include "lcio.h"



//using namespace lcio;



class VXDHelixFitterException : public std::exception {
   
   
protected:
   std::string message ;
   
   VXDHelixFitterException(){  /*no_op*/ ; } 
   
public: 
   virtual ~VXDHelixFitterException() { /*no_op*/; } 
   
   VXDHelixFitterException( const std::string& text ){
      message = "VXDHelixFitterException: " + text ;
   }
   
   virtual const char* what() const noexcept { return  message.c_str() ; } 
   
};



/** A class to make it quick to fit a track or hits and get back the chi2 and Ndf values and
 * also bundle the code used for that, so it doesn't have to be copied all over the places.
 * Uses a helix fit from the MarlinTrk class HelixFit.cc
 * It is named VXDHelixFitter, because it makes some assumptions about the hits, that come from them
 * being on the VXD. Specifically the errors passed to the helix fit are calculated on the assumption,
 * that du and dv are errors in the xy plane.
 * If this class is intended to be used for hits on different detectors, a careful redesign is necessary!
 */
class VXDHelixFitter{
   
   
public:
   
  VXDHelixFitter( edm4hep::Track* track ) ;
  VXDHelixFitter( std::vector < edm4hep::TrackerHit > trackerHits ) ;
   
   
   double getChi2(){ return _chi2; }
   int getNdf(){ return _Ndf; }
   
   float getOmega(){ return _omega; }
   float getTanLambda(){ return _tanLambda; }
   float getPhi0(){ return _phi0; }
   float getD0(){ return _d0; }
   float getZ0(){ return _z0; }
   
private:
   
  
   
   void fit();
   
   double _chi2;
   int _Ndf;
   
   float _omega;
   float _tanLambda;
   float _phi0;
   float _d0;
   float _z0;
   
   std::vector< edm4hep::TrackerHit > _trackerHits;
  
   
};

#endif
