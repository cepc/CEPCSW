#include <DataHelper/LCPlane3D.h>

#include <string>
#include <sstream>
#include <iostream>
#include <cmath>
#include <float.h>
#include <exception>
// #include <CLHEP/Vector/ThreeVector.h>

LCPlane3D::LCPlane3D(double a, double b, double c, double d) 
{
  _a = a;
  _b = b;
  _c = c;
  _d = d;
  normalize();

  // _d has to be negative to give the distance with the normal vector pointing
  // away from the Origen !!!
  if (_d > 0)
    {
      _a = -_a; 
      _b = -_b; 
      _c = -_c; 
      _d = -_d; 
    }
}

LCPlane3D::LCPlane3D(LCVector3D normal, LCVector3D point)
{
  LCVector3D n = normal.unit();

  _a = n.x();
  _b = n.y();
  _c = n.z();
  _d = - n*point;
}

LCPlane3D::LCPlane3D(LCVector3D point1, LCVector3D point2, LCVector3D point3) 
{
  LCVector3D n = ( (point2-point1).cross(point3-point1) ).unit();
  _a = n.x() ;
  _b = n.y() ;
  _c = n.z() ;

  _d = - n*point1;
}

LCPlane3D::LCPlane3D(LCVector3D normal, double distance) 
{
  LCVector3D n = normal.unit() ;
  _a = n.x();
  _b = n.y();
  _c = n.z();
  _d = -distance;
}

LCPlane3D::LCPlane3D(const LCPlane3D & plane) 
  : _a(plane._a),
    _b(plane._b),
    _c(plane._c),
    _d(plane._d)
{}

LCPlane3D & LCPlane3D::operator=(const LCPlane3D & rhs) 
{
  _a = rhs._a;
  _b = rhs._b;
  _c = rhs._c;
  _d = rhs._d;

  return *this;
}

double LCPlane3D::a() const 
{
  return _a ;
}

double LCPlane3D::b() const 
{
  return _b ;
}

double LCPlane3D::c() const 
{
  return _c ;
}

double LCPlane3D::d() const 
{
  return _d ;
}

LCVector3D LCPlane3D::normal() const 
{
  LCVector3D n(_a,_b,_c);
  return n.unit();
}

LCPlane3D & LCPlane3D::normalize() 
{
  double norm = sqrt( _a*_a + _b*_b + _c*_c );

  if (norm > 0.)
    {
      _a /= norm; 
      _b /= norm; 
      _c /= norm; 
      _d /= norm; 
    }

  return *this;
}

double LCPlane3D::distance(const LCVector3D & point) const 
{
  return _a*point.x() + _b*point.y() + _c*point.z() + _d ;
}

LCVector3D LCPlane3D::projectPoint(const LCVector3D & point) const 
{
  double k = distance(point) / ( _a*_a + _b*_b + _c*_c );
  return LCVector3D( point.x()-_a*k, point.y()-_b*k, point.z()-_c*k);
}

LCVector3D LCPlane3D::projectPoint() const 
{
  double k = -_d / ( _a*_a + _b*_b + _c*_c );
  return LCVector3D( _a*k, _b*k, _c*k);
}

bool LCPlane3D::operator==(const LCPlane3D & plane) const 
{
  return ( _a == plane._a && 
	   _b == plane._b && 
	   _c == plane._c && 
	   _d == plane._d );
}

bool LCPlane3D::operator!=(const LCPlane3D & plane) const 
{
  return ( _a != plane._a || 
	   _b != plane._b || 
	   _c != plane._c || 
	   _d != plane._d );
}

std::ostream & operator << (std::ostream &os, const LCPlane3D &p) 
{
  std::stringstream returnString ;
  bool isFirst = true;

  returnString << "(" ;
  if (p.a() != 0) 
    {
      returnString << p.a() << "*x" ;
      isFirst = false;
    }
  if (p.b() != 0)
    {
      if (!isFirst && (p.b() > 0.) ) returnString << "+";
      returnString << p.b() << "*y" ;
      isFirst = false;
    }
  if (p.c() != 0)
    {
      if (!isFirst && (p.c() > 0) ) returnString << "+";
      returnString << p.c() << "*z" ;
      isFirst = false;
    }
  if (p.d() != 0)
    {
      if (!isFirst && (p.d() > 0) ) returnString << "+";
      returnString << p.d() ;
    }
  returnString << "=0)" ;

  return os << returnString.str() ;
}
