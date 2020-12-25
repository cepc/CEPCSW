#ifndef INCLUDE_MarlinTrkUtils
#define INCLUDE_MarlinTrkUtils 1

#include <vector>
#include <array>
#include <cfloat>
#include <edm4hep/TrackerHitConst.h>

#include <LCIOSTLTypes.h>

namespace edm4hep{
  class Track;
  class TrackerHit;
  class TrackState;
}

namespace UTIL {
  class BitField64;
}


namespace MarlinTrk{
  class IMarlinTrack ;
}


namespace MarlinTrk{
  enum Location{
    AtOther = 0,
    AtIP = 1,
    AtFirstHit = 2,
    AtLastHit = 3,
    AtCalorimeter = 4,
    AtVertex = 5,
    LastLocation = 6
  };
  /** Takes a list of hits and uses the IMarlinTrack inferface to fit them using a supplied prefit containing
   *  a covariance matrix for the initialisation. The TrackImpl will have the 4 trackstates added to
   *  it @IP, @First_Hit, @Last_Hit and @CaloFace */
  int createFinalisedLCIOTrack(
      IMarlinTrack* marlinTrk,
      std::vector<edm4hep::ConstTrackerHit>& hit_list,
      edm4hep::Track* track,
      bool fit_backwards,
      edm4hep::TrackState* pre_fit,
      float bfield_z,
      double maxChi2Increment=DBL_MAX);
  
  /** Takes a list of hits and uses the IMarlinTrack inferface to fit them using a supplied covariance matrix
   *  for the initialisation. The TrackImpl will have the 4 trackstates added to
   *  it @IP, @First_Hit, @Last_Hit and @CaloFace */
  int createFinalisedLCIOTrack(
      IMarlinTrack* marlinTrk,
      std::vector<edm4hep::ConstTrackerHit>& hit_list,
      edm4hep::Track* track,
      bool fit_backwards,
      const std::array<float,15>& initial_cov_for_prefit,
      float bfield_z,
      double maxChi2Increment=DBL_MAX);
  
  /** Provides the values of a track state from the first, middle and last hits in the hit_list. */
  int createPrefit( std::vector<edm4hep::ConstTrackerHit>& hit_list, edm4hep::TrackState* pre_fit, float bfield_z, bool fit_backwards );

  /** Takes a list of hits and uses the IMarlinTrack inferface to fit them using a supplied prefit containing a covariance matrix for the initialisation. */  
  int createFit( std::vector<edm4hep::ConstTrackerHit>& hit_list, IMarlinTrack* marlinTrk, edm4hep::TrackState* pre_fit, float bfield_z, bool fit_backwards, double maxChi2Increment=DBL_MAX );

  /** Takes a fitted MarlinTrack, TrackImpl to record the fit and the hits which have been added to the fit.
   *  The TrackImpl will have the 4 trackstates added to it @IP, @First_Hit, @Last_Hit and @CaloFace.
   *  Note: the hit list is needed as the IMarlinTrack only contains the hits used in the fit, not the spacepoints
   *  (if any have been included) so as the strip hits cannot point to the space points we need to have the list so
   *  that they can be recorded in the LCIO TrackImpl */
  int finaliseLCIOTrack(
      IMarlinTrack* marlinTrk,
      edm4hep::Track* track,
      std::vector<edm4hep::ConstTrackerHit>& hit_list,
      edm4hep::TrackState* atLastHit=0,
      edm4hep::TrackState* atCaloFace=0);
  
  /** Set the subdetector hit numbers for the TrackImpl */
  void addHitNumbersToTrack(edm4hep::Track* track, std::vector<edm4hep::ConstTrackerHit>& hit_list, bool hits_in_fit, UTIL::BitField64& cellID_encoder);

  /** Set the subdetector hit numbers for the TrackImpl */
  void addHitNumbersToTrack(edm4hep::Track* track, std::vector<std::pair<edm4hep::ConstTrackerHit , double> >& hit_list, bool hits_in_fit, UTIL::BitField64& cellID_encoder);
  
}

#endif
