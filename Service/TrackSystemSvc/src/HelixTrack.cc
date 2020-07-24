
#include "TrackSystemSvc/IMarlinTrack.h"
#include "TrackSystemSvc/HelixTrack.h"
#include <cmath>
#include <TVector3.h>
#include <kaltest/THelicalTrack.h>
#include <edm4hep/Vector3d.h> //plcio/DoubleThree.h>
//#include "streamlog/streamlog.h"

// defines if s of the helix increases in the direction of x2 to x3 
bool HelixTrack::forwards = true;

HelixTrack::HelixTrack( const edm4hep::Vector3d& x1, const edm4hep::Vector3d& x2, const edm4hep::Vector3d& x3, double Bz, bool direction ){

  // Make a KalTest THelicalTrack
  TVector3 p1( x1[0], x1[1], x1[2] );
  TVector3 p2( x2[0], x2[1], x2[2] );
  TVector3 p3( x3[0], x3[1], x3[2] );
  /*
  std::cout << "debug: "  << "HelixTrack::HelixTrack Create from hits: \n " 
	    << "P1 x = " << p1.x() << " y = " << p1.y() << " z = " << p1.z() << " r = " << p1.Perp() << "\n " 
	    << "P2 x = " << p2.x() << " y = " << p2.y() << " z = " << p2.z() << " r = " << p2.Perp() << "\n " 
	    << "P3 x = " << p3.x() << " y = " << p3.y() << " z = " << p3.z() << " r = " << p3.Perp() << "\n "
	    << "Bz = " << Bz << " direction = " << direction 
	    << std::endl;
  */
  THelicalTrack*  helicalTrack;

  helicalTrack = new THelicalTrack( p1, p2, p3, Bz, direction ); 
    
  // Set the track parameters and convert from the KalTest system to the lcio system
  
  _phi0 = toBaseRange( helicalTrack->GetPhi0() + M_PI/2. ) ;
  _omega = 1. / helicalTrack->GetRho();
  _z0 = helicalTrack->GetDz();
  _d0 = - helicalTrack->GetDrho();
  _tanLambda = helicalTrack->GetTanLambda();
  
  _ref_point_x =  helicalTrack->GetPivot().X() ;
  _ref_point_y =  helicalTrack->GetPivot().Y() ;
  _ref_point_z =  helicalTrack->GetPivot().Z() ;
  
  delete helicalTrack;
  
}
  
  
HelixTrack::HelixTrack( const edm4hep::Vector3d& position, const edm4hep::Vector3d& p, double charge, double Bz ){
  
  _ref_point_x = position[0] ;
  _ref_point_y = position[1] ;
  _ref_point_z = position[2] ;
  
  _d0 = 0.0 ;
  _z0 = 0.0 ;
  
  const double pt = sqrt(p[0]*p[0]+p[1]*p[1]) ;
  
  double radius = pt / (2.99792458E-4*Bz) ; // for r in mm, p in GeV and Bz in Tesla
  
  _omega = charge/radius ;
  _tanLambda = p[2]/pt ;
  
  _phi0 = atan2(p[1],p[0]);
  
  _phi0 = toBaseRange(_phi0);

}

double HelixTrack::moveRefPoint( double x, double y, double z){
  
  const double radius = 1.0/_omega ; 
  
  const double sinPhi0 = sin(_phi0) ;
  const double cosPhi0 = cos(_phi0) ;
  
  const double deltaX = x - _ref_point_x ;
  const double deltaY = y - _ref_point_y ;
  
  double phi0Prime = atan2( sinPhi0 - (deltaX/(radius-_d0)) , cosPhi0 + (deltaY/(radius-_d0)) ) ;
  
  while ( phi0Prime < 0 )  phi0Prime += 2.0*M_PI ;
  while ( phi0Prime >= 2.0*M_PI ) phi0Prime -= 2.0*M_PI ;
  
  const double d0Prime = _d0 + deltaX*sinPhi0 - deltaY*cosPhi0 + ( ( deltaX*cosPhi0 + deltaY*sinPhi0 ) * tan( (phi0Prime-_phi0) / 2.0) ) ;
  
  // In order to have terms which behave well as Omega->0 we make use of deltaX and deltaY to replace sin( phi0Prime - phi0 ) and cos( phi0Prime - phi0 )
  
  const double sinDeltaPhi = ( -_omega / ( 1.0 - ( _omega * d0Prime ) ) ) * ( deltaX * cosPhi0 + deltaY * sinPhi0 ) ;
  
  const double cosDeltaPhi  = 1.0 + ( _omega*_omega / ( 2.0 * ( 1.0 - _omega * d0Prime ) ) ) * ( d0Prime*d0Prime - ( deltaX + _d0 * sinPhi0 )*( deltaX + _d0 * sinPhi0 ) - ( deltaY - _d0 * cosPhi0 )*( deltaY - _d0 * cosPhi0 ) ) ;
  
  const double s = atan2(-sinDeltaPhi,cosDeltaPhi) / _omega ;
  
  const double z0Prime  = _ref_point_z - z + _z0 + _tanLambda * s ;
  
  phi0Prime = toBaseRange(phi0Prime);

  _d0   = d0Prime ;
  _phi0 = phi0Prime ;
  _z0   = z0Prime ;

  _ref_point_x = x; 
  _ref_point_y = y;
  _ref_point_z = z;   
  
  return (s/radius);
  
}





