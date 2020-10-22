#include <DataHelper/LCLine3D.h>
#include <DataHelper/LCPlane3D.h>
#include <iostream>

#include <cmath>
#include <float.h>
#include <exception>

LCLine3D::LCLine3D()
{
  _reference.set(0.,0.,0.);
  _point.set(0.,0.,0.);
  _direction.set(1.,0.,0.);
}

LCLine3D::LCLine3D(const LCVector3D & point, const LCVector3D & direction) 
{
  set( point, direction, LCVector3D(0.,0.,0.) );
}

LCLine3D::LCLine3D(const LCVector3D & point,
		   const LCVector3D & direction,
		   const LCVector3D & reference) 
{
  set( point, direction, reference );
}

LCLine3D::LCLine3D(double d0, double phi0, double z0, double tanLambda) 
{
  set( d0, phi0, z0, tanLambda, LCVector3D(0.,0.,0.) );
}

LCLine3D::LCLine3D(double d0, double phi0, double z0, double tanLambda,
	 const LCVector3D & reference) 
{
  set( d0, phi0, z0, tanLambda, reference );
}

LCLine3D::LCLine3D(const LCLine3D & line) 
{
  _point     = line._point;
  _direction = line._direction;
  _reference = line._reference;
}

bool LCLine3D::set(const LCVector3D & point,
		   const LCVector3D & direction,
		   const LCVector3D & reference) 
{
  //  std::cout << "LCLine  " << _reference << " " << point << " " << direction << std::endl;

  _reference = reference;
  _direction = direction.unit();
  if (_direction.mag2() == 0)
    {
      return false;
    }

// calculate _point to be the PCA to the reference point according to the
// definition given in LC-LC-DET-2006-004:
// the x,y compnents have to beh teh  PCA, the z compnente is calculated 
// after that. 

  LCVector3D p = point, d = _direction;
  p.setZ(0.);
  d.setZ(0.); 

  if (d.mag() !=  0.)
    {
      double sFaktor = 1./d.mag();
      d = d.unit();
      double s = ( - p*d ) / d.mag2() ;
      // x,y componentes:
      //      _point = p + s * d ;
      // z component:
      _point = ( (point + s*sFaktor*_direction) );
    }
  else
    {
      _point = point;
      _point.setZ(0.);
    }
  //  std::cout << "LCLine: " << _reference << " " << _point << " " << _direction << std::endl;
  return true;
}

bool LCLine3D::set(double d0, double phi0, double z0, double tanLambda,
	 const LCVector3D & reference) 
{
  _reference = reference;
  _direction.set( cos(phi0), sin(phi0), tanLambda );
  _direction = _direction.unit();
  if (d0 == 0.)
    {
      _point.set(0.,0.,z0);
    }
  else
    {
      _point.set( ( d0*sin(phi0) ), ( d0*cos(phi0) ), z0 );
    }
  return true;
}

LCLine3D & LCLine3D::operator=(const LCLine3D & rhs) 
{
  _point     = rhs._point;
  _direction = rhs._direction;
  _reference = rhs._reference;

  return *this;
}

LCVector3D LCLine3D::position(const double s) const 
{
  return (_reference+_point + s*_direction) ;
}

LCVector3D LCLine3D::direction() const 
{
  return _direction;
}

double LCLine3D::distance(const LCVector3D & point) const 
{
  return ( point - position( projectPoint( point ) ) ).mag() ;
}

double LCLine3D::projectPoint(const LCVector3D & point) const 
{
  // the last therm : (...) / _direction.mag2() is not there becaus 
  // the _direction vector is normalised.
  //  return ( 2*point*_direction - _point*_direction ) / _direction.mag2() ;
  //  return ( point*_direction - (_reference+_point)*_direction ) / _direction.mag2() ;
  double x = ( point*_direction - (_reference+_point)*_direction ) / _direction.mag2() ;
  //  std::cout << "x: " << x << std::endl;
  //  std::cout << "point: " << point 
  //	    << " direction: " << _direction << " _point: " << _point 
  //	    << " d.mag: " << _direction.mag2() << std::endl;


  return x;
}

bool LCLine3D::operator==(const LCLine3D & rhs) const 
{
  return (_point == rhs._point &&
	  _direction == rhs._direction && 
	  _reference == rhs._reference);
}

bool LCLine3D::operator!=(const LCLine3D & rhs) const 
{
  return (_point != rhs._point ||
	  _direction != rhs._direction ||
	  _reference != rhs._reference) ;
}

double LCLine3D::intersectionWithPlane(const LCPlane3D plane, bool& pointExists) const 
{
  double c = direction() * plane.normal() ;

  if (c == 0)
    { // no interaction 
      pointExists = false;
      return DBL_MAX;
    }

  pointExists = true;
  return - ( position() * plane.normal() + plane.d() ) / c ;
}

std::ostream & operator << (std::ostream &os, const LCLine3D &l)
{
  return os << l.position() << "+s*" << l.direction() ;
}
