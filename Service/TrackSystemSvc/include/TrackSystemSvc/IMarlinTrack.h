#ifndef IMarlinTrack_h
#define IMarlinTrack_h

#include <cfloat>

//#include "lcio.h"

#include "edm4hep/TrackerHit.h"
#include "edm4hep/TrackerHitConst.h"
#include "edm4hep/TrackState.h"

//#include "gearimpl/Vector3D.h"
//#include "plcio/DoubleThree.h"
#include "edm4hep/Vector3d.h"

#include <exception>

  /** Interface for generic tracks in MarlinTrk. The interface should provide the functionality to
   *  perform track finding and fitting. It is asssumed that the underlying implemetation will by 
   *  a Kalman Filter or a similar algorithm.
   *
   * @version $Id: IMarlinTrack.h 3989 2012-09-04 07:16:44Z aplin $
   * @author S.Aplin, F. Gaede DESY
   */
namespace MarlinTrk{
  class IMarlinTrack {
    
  public:
    
    /** boolean constant for defining backward direction - to be used for intitialise */
    static const bool backward ;
    
    /** boolean constant for defining backward direction - to be used for intitialise */
    static const bool forward ;
    
    
    
    static const int modeBackward ;
    static const int modeClosest  ;
    static const int modeForward  ;
    
    
    static const int success ;  // no error
    static const int error ;
    static const int bad_intputs ;
    static const int no_intersection ; // no intersection found
    static const int site_discarded ;  // measurement discarded by the fitter
    static const int site_fails_chi2_cut ;  // measurement discarded by the fitter due to chi2 cut
    static const int all_sites_fail_fit ;   // no single measurement added to the fit
    
    
    /**default d'tor*/
    virtual ~IMarlinTrack() {};
    
    /** add hit to track - the hits have to be added ordered in time ( i.e. typically outgoing )
     *  this order will define the direction of the energy loss used in the fit
     */
    virtual int addHit(edm4hep::ConstTrackerHit hit) = 0 ;
    
    /** initialise the fit using the hits added up to this point -
     *  the fit direction has to be specified using IMarlinTrack::backward or IMarlinTrack::forward. 
     *  this is the order  wrt the order used in addHit() that will be used in the fit() 
     */
    virtual int initialise( bool fitDirection ) = 0 ; 
    
    /** initialise the fit with a track state, and z component of the B field in Tesla.
     *  the fit direction has to be specified using IMarlinTrack::backward or IMarlinTrack::forward. 
     *  this is the order that will be used in the fit().
     *  it is the users responsibility that the track state is consistent with the order
     *  of the hits used in addHit() ( i.e. the direction of energy loss )
     */
    virtual int initialise(  const edm4hep::TrackState& ts, double bfield_z, bool fitDirection ) = 0 ;
    
    
    /** perform the fit of all current hits, returns error code ( IMarlinTrack::success if no error ) .
     *  the fit will be performed  in the order specified at initialise() wrt the order used in addHit(), i.e.
     *  IMarlinTrack::backward implies fitting from the outside to the inside for tracks comming from the IP.
     */
    virtual int fit( double maxChi2Increment=DBL_MAX ) = 0 ;
    
    
    /** update the current fit using the supplied hit, return code via int. Provides the Chi2 increment to the fit from adding the hit via reference. 
     *  the given hit will not be added if chi2increment > maxChi2Increment. 
     */
    virtual int addAndFit( edm4hep::ConstTrackerHit hit, double& chi2increment, double maxChi2Increment=DBL_MAX ) = 0 ;

    
    /** obtain the chi2 increment which would result in adding the hit to the fit. This method will not alter the current fit, and the hit will not be stored in the list of hits or outliers
     */
    virtual int testChi2Increment( edm4hep::ConstTrackerHit hit, double& chi2increment ) = 0 ;

    
    /** smooth all track states 
     */
    virtual int smooth() = 0 ;
    
    
    /** smooth track states from the last filtered hit back to the measurement site associated with the given hit 
     */
    virtual int smooth( edm4hep::ConstTrackerHit hit ) = 0 ;
    
    
    
    // Track State Accessesors
    
    /** get track state, returning TrackState, chi2 and ndf via reference 
     */
    virtual int getTrackState( edm4hep::TrackState& ts, double& chi2, int& ndf ) = 0 ;
    
    
    /** get track state at measurement associated with the given hit, returning TrackState, chi2 and ndf via reference 
     */
    virtual int getTrackState( edm4hep::ConstTrackerHit hit, edm4hep::TrackState& ts, double& chi2, int& ndf ) = 0 ;
    
    /** get the list of hits included in the fit, together with the chi2 contributions of the hits. 
     *  Pointers to the hits together with their chi2 contribution will be filled into a vector of 
     *  pairs consitining of the pointer as the first part of the pair and the chi2 contribution as
     *  the second.
     */
    virtual int getHitsInFit( std::vector<std::pair<edm4hep::ConstTrackerHit, double> >& hits ) = 0 ;

    /** get the list of hits which have been rejected by from the fit due to the a chi2 increment greater than threshold,
     *  Pointers to the hits together with their chi2 contribution will be filled into a vector of 
     *  pairs consitining of the pointer as the first part of the pair and the chi2 contribution as
     *  the second.
     */
    virtual int getOutliers( std::vector<std::pair<edm4hep::ConstTrackerHit, double> >& hits ) = 0 ;
    
    /** get the current number of degrees of freedom for the fit.
     */
    virtual int getNDF( int& ndf ) = 0 ;
    
    /** get TrackeHit at which fit became constrained, i.e. ndf >= 0
     */
    virtual int getTrackerHitAtPositiveNDF( edm4hep::ConstTrackerHit& trkhit ) = 0 ;
    
    // PROPAGATORS 
    
    /** propagate the fit to the point of closest approach to the given point, returning TrackState, chi2 and ndf via reference    
     */
    virtual int propagate( const edm4hep::Vector3d& point, edm4hep::TrackState& ts, double& chi2, int& ndf ) = 0 ;
    
    
    /** propagate the fit at the measurement site associated with the given hit, to the point of closest approach to the given point,
     *  returning TrackState, chi2 and ndf via reference   
     */
    virtual int propagate( const edm4hep::Vector3d& point, edm4hep::ConstTrackerHit hit, edm4hep::TrackState& ts, double& chi2, int& ndf ) = 0 ;
    
    
    /** propagate fit to numbered sensitive layer, returning TrackState, chi2, ndf and integer ID of the intersected sensitive detector element via reference 
     */
    virtual int propagateToLayer( int layerID, edm4hep::TrackState& ts, double& chi2, int& ndf, int& detElementID, int mode=modeClosest ) = 0  ;
    
    /** propagate the fit at the measurement site associated with the given hit, to numbered sensitive layer, 
     *  returning TrackState, chi2, ndf and integer ID of the intersected sensitive detector element via reference 
     */
    virtual int propagateToLayer(int layerID, edm4hep::ConstTrackerHit hit, edm4hep::TrackState& ts, double& chi2, int& ndf, int& detElementID, int mode=modeClosest ) = 0;
    
    /** propagate the fit to sensitive detector element, returning TrackState, chi2 and ndf via reference
     */
    virtual int propagateToDetElement( int detElementID, edm4hep::TrackState& ts, double& chi2, int& ndf, int mode=modeClosest ) = 0  ;
    
    /** propagate the fit at the measurement site associated with the given hit, to sensitive detector element, 
     *  returning TrackState, chi2 and ndf via reference 
     */
    virtual int propagateToDetElement( int detEementID, edm4hep::ConstTrackerHit hit, edm4hep::TrackState& ts, double& chi2, int& ndf, int mode=modeClosest ) = 0  ;
    
    
    
    // EXTRAPOLATORS
    
    /** extrapolate the fit to the point of closest approach to the given point, returning TrackState, chi2 and ndf via reference   
     */
    virtual int extrapolate( const edm4hep::Vector3d& point, edm4hep::TrackState& ts, double& chi2, int& ndf ) = 0 ;
    
    /** extrapolate the fit at the measurement site associated with the given hit, to the point of closest approach to the given point, 
     *  returning TrackState, chi2 and ndf via reference   
     */
    virtual int extrapolate( const edm4hep::Vector3d& point, edm4hep::ConstTrackerHit hit, edm4hep::TrackState& ts, double& chi2, int& ndf ) = 0 ;
    
    /** extrapolate the fit to numbered sensitive layer, returning TrackState, chi2, ndf and integer ID of the intersected sensitive detector element via reference
     */
    virtual int extrapolateToLayer( int layerID, edm4hep::TrackState& ts, double& chi2, int& ndf, int& detElementID, int mode=modeClosest ) = 0  ;
    
    /** extrapolate the fit at the measurement site associated with the given hit, to numbered sensitive layer, 
     *  returning TrackState, chi2, ndf and integer ID of the intersected sensitive detector element via reference 
     */
    virtual int extrapolateToLayer( int layerID, edm4hep::ConstTrackerHit hit, edm4hep::TrackState& ts, double& chi2, int& ndf, int& detElementID, int mode=modeClosest ) = 0  ;
    
    /** extrapolate the fit to sensitive detector element, returning TrackState, chi2 and ndf via reference
     */
    virtual int extrapolateToDetElement( int detElementID, edm4hep::TrackState& ts, double& chi2, int& ndf, int mode=modeClosest ) = 0  ;
    
    /** extrapolate the fit at the measurement site associated with the given hit, to sensitive detector element, 
     *  returning TrackState, chi2 and ndf via reference 
     */
    virtual int extrapolateToDetElement( int detEementID, edm4hep::ConstTrackerHit hit, edm4hep::TrackState& ts, double& chi2, int& ndf, int mode=modeClosest ) = 0  ;
    
    
    // INTERSECTORS
    
    /** extrapolate the fit to numbered sensitive layer, returning intersection point in global coordinates and integer ID of the 
     *  intersected sensitive detector element via reference 
     */
    virtual int intersectionWithLayer( int layerID, edm4hep::Vector3d& point, int& detElementID, int mode=modeClosest ) = 0  ;
    
    /** extrapolate the fit at the measurement site associated with the given hit, to numbered sensitive layer,
     *  returning intersection point in global coordinates and integer ID of the intersected sensitive detector element via reference 
     */
    virtual int intersectionWithLayer( int layerID, edm4hep::ConstTrackerHit hit, edm4hep::Vector3d& point, int& detElementID, int mode=modeClosest ) = 0  ;
    
    
    /** extrapolate the fit to numbered sensitive detector element, returning intersection point in global coordinates via reference 
     */
    virtual int intersectionWithDetElement( int detElementID, edm4hep::Vector3d& point, int mode=modeClosest ) = 0  ;
    
    /** extrapolate the fit at the measurement site associated with the given hit, to sensitive detector element,
     *  returning intersection point in global coordinates via reference 
     */
    virtual int intersectionWithDetElement( int detEementID, edm4hep::ConstTrackerHit hit, edm4hep::Vector3d& point, int mode=modeClosest ) = 0  ;
    
    
  protected:
    
  private:
    
    IMarlinTrack& operator=( const IMarlinTrack&) ; // disallow assignment operator 
                
  } ;
  std::string errorCode( int error );
}
#endif

