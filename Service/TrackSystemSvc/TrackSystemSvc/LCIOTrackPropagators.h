#ifndef LCIOTrackPropagators_h
#define LCIOTrackPropagators_h

#include <edm4hep/TrackState.h>

namespace LCIOTrackPropagators{
  
  /** Propagate trackstate to a new reference point
   */
  int PropagateLCIOToNewRef( edm4hep::TrackState& ts, double xref, double yref, double zref) ;
  
  /** Propagate trackstate to a new reference point taken as its crossing point with a cylinder of infinite length centered at x0,y0, parallel to the z axis. 
   For direction== 0  the closest crossing point will be taken
   For direction== 1  the first crossing traversing in positive s will be taken
   For direction==-1  the first crossing traversing in negative s will be taken
   */
  int PropagateLCIOToCylinder( edm4hep::TrackState& ts, float r, float x0, float y0, int direction=0, double epsilon=1.0e-8) ;
  
  
  /** Propagate trackstate to a new reference point taken as its crossing point with an infinite plane located at z, perpendicular to the z axis 
   */
  int PropagateLCIOToZPlane( edm4hep::TrackState& ts, float z) ;
  
  /** Propagate trackstate to a new reference point taken as its crossing point with a plane parallel to the z axis, containing points x1,x2 and y1,y2. Tolerance for intersection determined by epsilon.
   For direction== 0  the closest crossing point will be taken
   For direction== 1  the first crossing traversing in positive s will be taken
   For direction==-1  the first crossing traversing in negative s will be taken
   */
  int PropagateLCIOToPlaneParralelToZ( edm4hep::TrackState& ts, float x1, float y1, float x2, float y2, int direction=0, double epsilon=1.0e-8) ;
  
  
  
}

#endif
