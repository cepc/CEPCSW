#ifndef HelixTrack_h
#define HelixTrack_h

#include <cmath>

namespace edm4hep{
  class Vector3d;
}

class HelixTrack {
  
public:
  
  HelixTrack( double ref_point_x, double ref_point_y, double ref_point_z, double d0, double z0, double phi0, double omega, double tanLambda )
  : _ref_point_x(ref_point_x), _ref_point_y(ref_point_y), _ref_point_z(ref_point_z), _d0(d0), _z0(z0), _phi0(phi0), _omega(omega), _tanLambda(tanLambda)
  {    
    while ( _phi0 < -M_PI ) _phi0 += 2.0*M_PI ;
    while ( _phi0 >= M_PI ) _phi0 -= 2.0*M_PI;
  } 
  
  HelixTrack( const edm4hep::Vector3d& x1, const edm4hep::Vector3d& x2, const edm4hep::Vector3d& x3, double Bz, bool direction );
  
  HelixTrack( const edm4hep::Vector3d& position, const edm4hep::Vector3d& p, double charge, double Bz ) ;
  
  double moveRefPoint( double x, double y, double z) ;
  
  double  getRefPointX() const { return _ref_point_x ; }        
  double  getRefPointY() const { return _ref_point_y ; }        
  double  getRefPointZ() const { return _ref_point_z ; }        
  double  getD0() const { return _d0 ; }        
  double  getZ0() const { return _z0 ; }        
  double  getPhi0() const { return  _phi0; }    
  double  getOmega() const { return _omega ; }  
  double  getTanLambda() const { return _tanLambda ; } 
  
  // defines if s of the helix increases in the direction of x2 to x3 
  static bool forwards;
  
private:
  
  double _ref_point_x ;
  double _ref_point_y ;
  double _ref_point_z ;
  double _d0 ;
  double _z0 ;
  double _phi0 ;
  double _omega ;
  double _tanLambda ;
  
  /** helper function to restrict the range of the azimuthal angle to ]-pi,pi]*/
  inline double toBaseRange( double phi) const {
    while( phi <= -M_PI ){  phi += 2. * M_PI ; }
    while( phi >   M_PI ){  phi -= 2. * M_PI ; }
    return phi ;
  }
  
} ;




#endif
