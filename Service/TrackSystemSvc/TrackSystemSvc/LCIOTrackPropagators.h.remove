#ifndef LCIOTrackPropagators_h
#define LCIOTrackPropagators_h

namespace EVENT{
  class TrackState ;
}

namespace IMPL{
  class TrackStateImpl ;
}


namespace LCIOTrackPropagators{
  
  /** Propagate trackstate to a new reference point
   */
  int PropagateLCIOToNewRef( IMPL::TrackStateImpl& ts, double xref, double yref, double zref) ;
  
  /** Propagate trackstate to a new reference point taken as its crossing point with a cylinder of infinite length centered at x0,y0, parallel to the z axis. 
   For direction== 0  the closest crossing point will be taken
   For direction== 1  the first crossing traversing in positive s will be taken
   For direction==-1  the first crossing traversing in negative s will be taken
   */
  int PropagateLCIOToCylinder( IMPL::TrackStateImpl& ts, float r, float x0, float y0, int direction=0, double epsilon=1.0e-8) ;
  
  
  /** Propagate trackstate to a new reference point taken as its crossing point with an infinite plane located at z, perpendicular to the z axis 
   */
  int PropagateLCIOToZPlane( IMPL::TrackStateImpl& ts, float z) ;
  
  /** Propagate trackstate to a new reference point taken as its crossing point with a plane parallel to the z axis, containing points x1,x2 and y1,y2. Tolerance for intersection determined by epsilon.
   For direction== 0  the closest crossing point will be taken
   For direction== 1  the first crossing traversing in positive s will be taken
   For direction==-1  the first crossing traversing in negative s will be taken
   */
  int PropagateLCIOToPlaneParralelToZ( IMPL::TrackStateImpl& ts, float x1, float y1, float x2, float y2, int direction=0, double epsilon=1.0e-8) ;
  
  
  
}

#endif
