
#include "TrackSystemSvc/LCIOTrackPropagators.h"

#include <cmath>
#include <iostream>

#include "edm4hep/TrackState.h"

#include "CLHEP/Matrix/Matrix.h"
#include "CLHEP/Matrix/SymMatrix.h"
#include "CLHEP/Vector/TwoVector.h"

//#include "streamlog/streamlog.h"


namespace LCIOTrackPropagators{
  
  int PropagateLCIOToNewRef( edm4hep::TrackState& ts, double xref, double yref, double zref ) {
    
    //    std::cout << "PropagateLCIOToNewRef: x:y:z = " << xref << " : " << yref << " : " << zref << std::endl ;
    
    // Convert Parameters
    
    const double d0    = ts.D0 ;
    const double phi0  = ts.phi ;
    const double omega = ts.omega ;
    const double z0    = ts.Z0 ;
    const double tanL  = ts.tanLambda ;
    
    //   const double charge = omega/fabs(omega) ;
    edm4hep::Vector3f ref = ts.referencePoint ;  
    
    const double radius = 1.0/omega ; 
    
    const double sinPhi0 = sin(phi0) ;
    const double cosPhi0 = cos(phi0) ;
    
    const double deltaX = xref - ref[0] ;
    const double deltaY = yref - ref[1] ;
    
    double phi0Prime = atan2( sinPhi0 - (deltaX/(radius-d0)) , cosPhi0 + (deltaY/(radius-d0)) ) ;
    
    while ( phi0Prime < 0 )  phi0Prime += 2.0*M_PI ;
    while ( phi0Prime >= 2.0*M_PI ) phi0Prime -= 2.0*M_PI ;
    
    const double d0Prime = d0 + deltaX*sinPhi0 - deltaY*cosPhi0 + ( ( deltaX*cosPhi0 + deltaY*sinPhi0 ) * tan( (phi0Prime-phi0) / 2.0) ) ;
    
    // In order to have terms which behave well as Omega->0 we make use of deltaX and deltaY to replace sin( phi0Prime - phi0 ) and cos( phi0Prime - phi0 )
    
    const double sinDeltaPhi = ( -omega / ( 1.0 - ( omega * d0Prime ) ) ) * ( deltaX * cosPhi0 + deltaY * sinPhi0 ) ;
    
    const double cosDeltaPhi  = 1.0 + ( omega*omega / ( 2.0 * ( 1.0 - omega * d0Prime ) ) ) * ( d0Prime*d0Prime - ( deltaX + d0 * sinPhi0 )*( deltaX + d0 * sinPhi0 ) - ( deltaY - d0 * cosPhi0 )*( deltaY - d0 * cosPhi0 ) ) ;
    
    const double s = atan2(-sinDeltaPhi,cosDeltaPhi) / omega ;
    
    const double z0Prime  = ref[2] - zref + z0 + tanL * s ;
    
    
    // Convert Covariance Matrix
    CLHEP::HepSymMatrix cov0(5) ; 
    
    int icov = 0 ;
    
    for(int irow=0; irow<5; ++irow ){
      for(int jcol=0; jcol<irow+1; ++jcol){
        //      std::cout << "row = " << irow << " col = " << jcol << std::endl ;
        //      std::cout << "cov["<< icov << "] = " << _cov[icov] << std::endl ;
        cov0[irow][jcol] = ts.covMatrix[icov] ;
        ++icov ;
      }
    }
    
    CLHEP::HepMatrix propagatorMatrix(5, 5, 0) ;
    
    // LC_0 = { d0, phi0, omega, z0, tanLambda }
    
    // d d0' / d LC_0
    propagatorMatrix(1,1) = cosDeltaPhi ;
    propagatorMatrix(1,2) = -( radius - d0 ) * sinDeltaPhi ;
    propagatorMatrix(1,3) = radius*radius * ( cosDeltaPhi -1 ) ;
    
    // d phi0' / d LC_0 
    propagatorMatrix(2,1) = sinDeltaPhi / ( radius - d0Prime ) ;
    propagatorMatrix(2,2) = ( ( radius - d0 ) * cosDeltaPhi ) / ( radius - d0Prime ) ;
    propagatorMatrix(2,3) = radius*radius * sinDeltaPhi / ( radius - d0Prime ) ;
    
    // d omega' / d LC_0 
    propagatorMatrix(3,3) = 1.0 ;
    
    // d z0' / d LC_0 
    propagatorMatrix(4,1) = radius * tanL * sinDeltaPhi / ( d0Prime + radius ) ;
    propagatorMatrix(4,2) = radius * tanL * ( 1.0 - ( ( d0 + radius ) * cosDeltaPhi / ( d0Prime + radius ) ) ) ;
    propagatorMatrix(4,3) = radius*radius * tanL * ( (phi0Prime - phi0) - radius * sinDeltaPhi / ( d0Prime + radius ) ) ;
    propagatorMatrix(4,4) = 1.0 ;
    propagatorMatrix(4,5) = s ;
    
    // d tanLambda' / d LC_0 
    propagatorMatrix(5,5) = 1.0 ;
    
    
    CLHEP::HepSymMatrix covPrime =  cov0.similarity(propagatorMatrix);
    
    std::array<float,15> cov;
    
    icov = 0 ;
    
    for(int irow=0; irow<5; ++irow ){
      for(int jcol=0; jcol<irow+1; ++jcol){
        //      std::cout << "row = " << irow << " col = " << jcol << std::endl ;
        cov[icov] = covPrime[irow][jcol] ;
        //      std::cout << "lcCov["<< icov << "] = " << lcCov[icov] << std::endl ;
        ++icov ;
      }
    }
    
    while ( phi0Prime < -M_PI  ) phi0Prime += 2.0*M_PI ;
    while ( phi0Prime >= M_PI )  phi0Prime -= 2.0*M_PI ;
    
    ts.D0 = d0Prime;  
    ts.phi = phi0Prime ; 
    ts.omega = omega;
    ts.Z0 = z0Prime ;  
    ts.tanLambda = tanL;  
    
    
    float refPointPrime[3] ;
    refPointPrime[0] = xref ;
    refPointPrime[1] = yref ;
    refPointPrime[2] = zref ;
    
    ts.referencePoint = refPointPrime;
    
    ts.covMatrix = cov;
    
    return 0 ;
    
    
  }
  
  // Propagate track to a new reference point taken as its crossing point with a cylinder of infinite length centered at x0,y0, parallel to the z axis. 
  // For direction== 0  the closest crossing point will be taken
  // For direction== 1  the first crossing traversing in positive s will be taken
  // For direction==-1  the first crossing traversing in negative s will be taken
  
  int PropagateLCIOToCylinder( edm4hep::TrackState& ts, float r0, float x0, float y0, int direction, double epsilon){
    
    // taken from http://paulbourke.net/geometry/2circle/tvoght.c
    
    //    std::cout << "PropagateLCIOToCylinder: r = " << r0 << " x0:y0 = " << x0 << " : " << y0 << " direction = " << direction << std::endl ;
    
    
    const double x_ref = ts.referencePoint[0] ; 
    const double y_ref = ts.referencePoint[1] ; 
    const double z_ref = ts.referencePoint[2] ; 
    
    const double d0    = ts.D0 ;
    const double z0    = ts.Z0 ;
    const double phi0  = ts.phi ;
    const double tanl  = ts.tanLambda ;
    const double omega = ts.omega ;
    
    const double rho   = 1.0 / omega ;
    const double x_pca = x_ref - d0 * sin(phi0) ; 
    const double y_pca = y_ref + d0 * cos(phi0) ; 
    const double z_pca = z_ref + z0 ;
    
    const double sin_phi0 = sin(phi0); 
    const double cos_phi0 = cos(phi0); 
    
    const double x_c   = x_ref + ( rho - d0) * sin_phi0 ;
    const double y_c   = y_ref - ( rho - d0) * cos_phi0 ;
    
    const double r1 = fabs(rho) ;
    
    /* dx and dy are the vertical and horizontal distances between
     * the circle centers.
     */
    const double dx = x_c - x0;
    const double dy = y_c - y0;
    
    /* Determine the straight-line distance between the centers. */
    const double  d = hypot(dx,dy); 
    
    /* Check for solvability. */
    if (d > (r0 + r1))
      {
      /* no solution. circles do not intersect. */
      return 1;
      }
    if (d < fabs(r0 - r1))
      {
      /* no solution. one circle is contained in the other */
      return 2;
      }
    if (d < epsilon)
      {
      /* no solution. circles have common centre */
      return 3;
      }
    
    /* 'point 2' is the point where the line through the circle
     * intersection points crosses the line between the circle
     * centers.  
     */
    
    /* Determine the distance from point 0 to point 2. */
    const double a = ((r0*r0) - (r1*r1) + (d*d)) / (2.0 * d) ;
    
    /* Determine the coordinates of point 2. */
    const double x2 = x0 + (dx * a/d);
    const double y2 = y0 + (dy * a/d);
    
    /* Determine the distance from point 2 to either of the
     * intersection points.
     */
    const double h = sqrt((r0*r0) - (a*a));
    
    /* Now determine the offsets of the intersection points from
     * point 2.
     */
    const double rx = -dy * (h/d);
    const double ry =  dx * (h/d);
    
    /* Determine the absolute intersection points. */
    const double x_ins1 = x2 + rx;
    const double y_ins1 = y2 + ry;
    
    const double x_ins2 = x2 - rx;
    const double y_ins2 = y2 - ry;
    
    //    std::cout << "PropagateLCIOToCylinder: 1st solution x:y = " << x_ins1 << " : " << y_ins1 << std::endl ;
    //    std::cout << "PropagateLCIOToCylinder: 2nd solution x:y = " << x_ins2 << " : " << y_ins2 << std::endl ;
    
    // now calculate the path lengths 
    double s_1 = 0.0 ;
    double s_2 = 0.0 ;
    
    const double delta_x1 = x_ins1 - x_pca ;
    const double delta_y1 = y_ins1 - y_pca ;
    
    const double sin_delta_phi1 =     - omega*delta_x1*cos_phi0 - omega*delta_y1*sin_phi0 ;
    const double cos_delta_phi1 = 1.0 - omega*delta_x1*sin_phi0 + omega*delta_y1*cos_phi0 ;
    
    s_1 = atan2(-sin_delta_phi1,cos_delta_phi1) / omega ;
    
    
    const double delta_x2 = x_ins2 - x_pca ;
    const double delta_y2 = y_ins2 - y_pca ;
    
    const double sin_delta_phi2 =     - omega*delta_x2*cos_phi0 - omega*delta_y2*sin_phi0 ;
    const double cos_delta_phi2 = 1.0 - omega*delta_x2*sin_phi0 + omega*delta_y2*cos_phi0 ;
    
    s_2 = atan2(-sin_delta_phi2,cos_delta_phi2) / omega ;
    
    double x=0, y=0, z=0 ;
    
    if( direction == 0 ) { // take closest intersection
      if( fabs(s_1) < fabs(s_2) ) {
        //      std::cout << "PropagateLCIOToCylinder: use 1st solution" << std::endl ;
        x = x_ins1 ;
        y = y_ins1 ;
        z = z_pca + s_1 * tanl ;
      }
      else {
        //      std::cout << "PropagateLCIOToCylinder: use 2nd solution" << std::endl ;
        x = x_ins2 ;
        y = y_ins2 ;
        z = z_pca + s_2 * tanl ;          
      }
    }
    
    else {
      
      if ( s_1 < 0.0 ) s_1 +=  2.0*M_PI * r1 ;
      if ( s_2 < 0.0 ) s_2 +=  2.0*M_PI * r1 ;
      
      if( direction == 1 ){ // take the intersection with smallest s
        if( s_1 < s_2 ) {
          x = x_ins1 ;
          y = y_ins1 ;
          z = z_pca + s_1 * tanl ;
        }
        else {
          x = x_ins2 ;
          y = y_ins2 ;
          z = z_pca + s_2 * tanl ;        
        }
      } 
      else if(direction == -1) {  // else take the intersection with largest s 
        if( s_1 > s_2 ){
          x = x_ins1 ;
          y = y_ins1 ;
          z = z_pca + s_1 * tanl ;
        }
        else{
          x = x_ins2 ;
          y = y_ins2 ;
          z = z_pca + s_2 * tanl ;
        }
      }
    }
    
    
    return PropagateLCIOToNewRef(ts,x,y,z);
    
  }
  
  int PropagateLCIOToZPlane( edm4hep::TrackState& ts, float z) {
    
    
    const double x_ref = ts.referencePoint[0] ; 
    const double y_ref = ts.referencePoint[1] ; 
    const double z_ref = ts.referencePoint[2] ; 
    
    const double d0    = ts.D0 ;
    const double z0    = ts.Z0 ;
    const double phi0  = ts.phi ;
    const double tanl  = ts.tanLambda ;
    const double omega = ts.omega ;
    
    const double x_pca = x_ref - d0 * sin(phi0) ; 
    const double y_pca = y_ref + d0 * cos(phi0) ; 
    const double z_pca = z_ref + z0 ;
    
    // get path length to crossing point 
    const double s = ( z - z_pca ) / tanl;
    
    const double delta_phi_half = (omega*s)/2.0 ;
    
    const double x = x_pca + s * ( sin(delta_phi_half) / delta_phi_half ) *  cos( phi0 - delta_phi_half ) ;
    const double y = y_pca + s * ( sin(delta_phi_half) / delta_phi_half ) *  sin( phi0 - delta_phi_half ) ;
    
    return PropagateLCIOToNewRef(ts,x,y,z);
    
  }
  
  
  
  // Propagate track to a new reference point taken as its crossing point with a plane parallel to the z axis, containing points x1,x2 and y1,y2. Tolerance for intersection determined by epsilon.
  // For direction ==  0  the closest crossing point will be taken
  // For direction ==  1  the first crossing traversing in positive s will be taken
  // For direction == -1  the first crossing traversing in negative s will be taken
  int PropagateLCIOToPlaneParralelToZ( edm4hep::TrackState& ts, float x1, float y1, float x2, float y2, int direction, double epsilon) {
    
    // check that direction has one of the correct values
    if( !( direction == 0 || direction == 1 || direction == -1) ) return -1 ;
    
    // taken from http://paulbourke.net/geometry/sphereline/raysphere.c
    
    const double x_ref = ts.referencePoint[0] ; 
    const double y_ref = ts.referencePoint[1] ; 
    const double z_ref = ts.referencePoint[2] ; 
    
    const double d0    = ts.D0 ;
    const double z0    = ts.Z0 ;
    const double phi0  = ts.phi ;
    const double tanl  = ts.tanLambda ;
    const double omega = ts.omega ;
    
    const double rho   = 1.0 / omega ;
    const double x_pca = x_ref - d0 * sin(phi0) ; 
    const double y_pca = y_ref + d0 * cos(phi0) ; 
    const double z_pca = z_ref + z0 ;
    
    const double sin_phi0 = sin(phi0); 
    const double cos_phi0 = cos(phi0); 
    
    const double x_c   = x_ref + ( rho - d0) * sin_phi0 ;
    const double y_c   = y_ref - ( rho - d0) * cos_phi0 ;
    
    const double dx = x2 - x1 ;
    const double dy = y2 - y1 ;
    
    const double a = dx * dx + dy * dy ;
    
    const double b = 2.0 * ( dx * (x1 - x_c) + dy * (y1 - y_c) ) ;
    
    double c = x_c * x_c + y_c * y_c;
    c += x1 * x1 + y1 * y1 ;
    
    c -= 2.0 * ( x_c * x1 + y_c * y1 );
    c -= rho * rho;
    
    const double bb4ac = b * b - 4.0 * a * c;
    
    double u1 ; // first solution
    double u2 ; // second solution
    
    /* Check for solvability. */
    if (bb4ac + epsilon < 0.0 ) { // no intersection
      return(1);
    }
    else if(bb4ac - epsilon < 0.0) { // circle intersects at one point, tangential 
      return(2);
    }
    else{
      u1 = (-b + sqrt(bb4ac)) / (2.0 * a);
      u2 = (-b - sqrt(bb4ac)) / (2.0 * a);
    }
    
    const double x_ins1 = x1 + u1 * (dx) ;
    const double y_ins1 = y1 + u1 * (dy) ;
    
    const double x_ins2 = x1 + u2 * (dx) ;
    const double y_ins2 = y1 + u2 * (dy) ;
    
    //    std::cout << "PropagateLCIOToPlaneParralelToZ: 1st solution u = " << u1 << " : " << u1 * dx << " " << u1 * dy << std::endl ;  
    //    std::cout << "PropagateLCIOToPlaneParralelToZ: 2nd solution u = " << u2 << " : " << u2 * dx << " " << u2 * dy << std::endl ;
    //    std::cout << "PropagateLCIOToPlaneParralelToZ: 1st solution x:y = " << x_ins1 << " : " << y_ins1 << std::endl ;
    //    std::cout << "PropagateLCIOToPlaneParralelToZ: 2nd solution x:y = " << x_ins2 << " : " << y_ins2 << std::endl ;
    
    // now calculate the path lengths 
    double s_1 = 0.0 ;
    double s_2 = 0.0 ;
    
    const double delta_x1 = x_ins1 - x_pca ;
    const double delta_y1 = y_ins1 - y_pca ;
    
    const double sin_delta_phi1 =     - omega*delta_x1*cos_phi0 - omega*delta_y1*sin_phi0 ;
    const double cos_delta_phi1 = 1.0 - omega*delta_x1*sin_phi0 + omega*delta_y1*cos_phi0 ;
    
    s_1 = atan2(-sin_delta_phi1,cos_delta_phi1) / omega ;
    
    const double delta_x2 = x_ins2 - x_pca ;
    const double delta_y2 = y_ins2 - y_pca ;
    
    const double sin_delta_phi2 =     - omega*delta_x2*cos_phi0 - omega*delta_y2*sin_phi0 ;
    const double cos_delta_phi2 = 1.0 - omega*delta_x2*sin_phi0 + omega*delta_y2*cos_phi0 ;
    
    s_2 = atan2(-sin_delta_phi2,cos_delta_phi2) / omega ;
    
    double x(0.0), y(0.0), z(0.0) ;
    
    if( direction == 0 ) { // take closest intersection
      if( fabs(s_1) < fabs(s_2) ) {
        //      std::cout << "PropagateLCIOToPlaneParralelToZ: take 1st solution " << std::endl;
        x = x_ins1 ;
        y = y_ins1 ;
        z = z_pca + s_1 * tanl ;
      }
      else {
        //      std::cout << "PropagateLCIOToPlaneParralelToZ: take 2nd solution " << std::endl;
        x = x_ins2 ;
        y = y_ins2 ;
        z = z_pca + s_2 * tanl ;          
      }
    }
    
    else{
      
      if ( s_1 < 0.0 ) s_1 +=  2.0*M_PI * fabs(rho) ;
      if ( s_2 < 0.0 ) s_2 +=  2.0*M_PI * fabs(rho) ;
      
      if( direction == 1 ){ // take the intersection with smallest s
        if( s_1 < s_2 ) {
          //      std::cout << "PropagateLCIOToPlaneParralelToZ: take 1st solution " << std::endl;
          x = x_ins1 ;
          y = y_ins1 ;
          z = z_pca + s_1 * tanl ;
        }
        else {
          //      std::cout << "PropagateLCIOToPlaneParralelToZ: take 2nd solution " << std::endl;
          x = x_ins2 ;
          y = y_ins2 ;
          z = z_pca + s_2 * tanl ;        
        }
      } 
      else if( direction == -1 ) {  // else take the intersection with largest s 
        if( s_1 > s_2 ){
          //      std::cout << "PropagateLCIOToPlaneParralelToZ: take 1st solution " << std::endl;
          x = x_ins1 ;
          y = y_ins1 ;
          z = z_pca + s_1 * tanl ;
        }
        else{
          //      std::cout << "PropagateLCIOToPlaneParralelToZ: take 2nd solution " << std::endl;
          x = x_ins2 ;
          y = y_ins2 ;
          z = z_pca + s_2 * tanl ;
        }
      }
    }
    
    return PropagateLCIOToNewRef(ts,x,y,z);
    
  }
  
  
} // end of TrackPropagators
