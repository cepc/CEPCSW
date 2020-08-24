
#include "TrackSystemSvc/IMarlinTrack.h"

namespace MarlinTrk{

  const bool IMarlinTrack::backward = false ;
  const bool IMarlinTrack::forward  = ! IMarlinTrack::backward;
  
  const int IMarlinTrack::modeBackward = - 1 ;
  const int IMarlinTrack::modeClosest  =   0 ;
  const int IMarlinTrack::modeForward  = + 1 ;
  
  
  const int IMarlinTrack::success  = 0 ;  // no error
  const int IMarlinTrack::error = 1 ;
  const int IMarlinTrack::bad_intputs = 3 ;
  const int IMarlinTrack::no_intersection = 4 ; // no intersection found
  const int IMarlinTrack::site_discarded = 5 ;  // measurement discarded by the fitter
  const int IMarlinTrack::site_fails_chi2_cut = 6 ;  // measurement discarded by the fitter due to chi2 cut
  const int IMarlinTrack::all_sites_fail_fit = 7 ;   // no single measurement added to the fit
  
  
  /** Helper function to convert error return code to string */
  std::string errorCode( int error ){
    
    switch ( error ){
      case IMarlinTrack::success                : return "IMarlinTrack::success";             break;
      case IMarlinTrack::error                  : return "IMarlinTrack::error";               break;
      case IMarlinTrack::bad_intputs            : return "IMarlinTrack::bad_intputs";         break;
      case IMarlinTrack::no_intersection        : return "IMarlinTrack::no_intersection";     break;
      case IMarlinTrack::site_discarded         : return "IMarlinTrack::site_discarded";      break;
      case IMarlinTrack::site_fails_chi2_cut    : return "IMarlinTrack::site_fails_chi2_cut"; break;
      case IMarlinTrack::all_sites_fail_fit     : return "IMarlinTrack::all_sites_fail_fit";  break;
      default: return "UNKNOWN" ;
    }
  }
  
}
