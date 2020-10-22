#include "DataHelper/SimpleHelix.h"
#include "CLHEP/Matrix/Matrix.h"
#include "CLHEP/Matrix/Vector.h"

#include <iostream>
#include <iomanip>
#include <DataHelper/LCLine3D.h>
#include <cmath>
#include <float.h>
#include <exception>

const double SimpleHelix::_a  = 2.99792458E-4;
const double SimpleHelix::_pi = M_PI;

SimpleHelix::SimpleHelix( double d0, double phi0, double omega,
			  double z0, double tanLambda,
			  LCVector3D referencePoint, LCErrorMatrix* errors) 
{
  init();

  _d0        = d0;
  _phi0      = phi0;
  _omega     = omega;
  _z0        = z0;
  _tanLambda = tanLambda;

  _reference = referencePoint;

  //  _Bz = Bz;

  if ( errors == 0 )
    {
      _errors = new LCErrorMatrix(5,0);
    }
  else
    {
      if ( errors->num_row() != 5 )
	throw("The error matrix has to be a 5x5 symmetric matrix in SimpleHelix Constructor");
      *_errors = *errors;
    }
}

void SimpleHelix::init()
{
  _d0         = 0;
  _phi0       = 0;
  _omega      = 1;
  _z0         = 0;
  _tanLambda  = 0;
  //  _Bz         = 4;

  _reference.set(0.,0.,0.); 

  _helixStart = -DBL_MAX;
  _helixEnd   =  DBL_MAX;

}

double SimpleHelix::getCentreX() const
{
  return ( _reference.x() + ((1/_omega) - _d0) * sin(_phi0) );
}

double SimpleHelix::getCentreY() const
{
  return ( _reference.y() - ((1/_omega) - _d0) * cos(_phi0) );
}

double SimpleHelix::getWindingLength() const
{
  return 2*_pi*sqrt(1+_tanLambda*_tanLambda)/fabs(_omega);
}

LCVector3D SimpleHelix::getPosition(double s, LCErrorMatrix* errors) const 
{
  LCVector3D x;

  double xc = getCentreX();
  double yc = getCentreY();

  double varphi0 = _phi0 + ((_omega * _pi) / (2*fabs(_omega)));
  double w =  _omega / sqrt(1 + _tanLambda*_tanLambda);

  x.setX(xc + fabs(1/_omega) * cos( w*s - varphi0 ) );
  x.setY(yc + fabs(1/_omega) * sin( (-w*s) + varphi0) );
  x.setZ(_reference.z() + _z0 + s*_tanLambda/sqrt(1 + _tanLambda*_tanLambda) );

  return x;
}

LCVector3D SimpleHelix::getDirection(double s,  LCErrorMatrix* errors) const
{
  LCVector3D t;

  double varphi0 = _phi0 + ((_omega * _pi) / (2*fabs(_omega)));
  double w =  _omega / sqrt(1 + _tanLambda*_tanLambda);

  t.setX(-w*fabs(1/_omega) * sin( (w*s) - varphi0) );
  t.setY(-w*fabs(1/_omega) * cos( (-w*s) + varphi0) );
  t.setZ(_tanLambda/sqrt(1 + _tanLambda*_tanLambda) );

  return t.unit();
}

LCErrorMatrix SimpleHelix::getCovarianceMatrix( double s) const  
{
  return LCErrorMatrix( 6 , 0 ) ;
}

double SimpleHelix::getPathAt(const LCVector3D position ) const
{
  double sStart = (position.z() - _reference.z() - _z0)
    * sqrt(1 + _tanLambda*_tanLambda) / _tanLambda;

  if (sStart <= _helixStart) sStart = _helixStart;
  if (sStart >= _helixEnd) sStart = _helixEnd;

  double startRange = sStart - getWindingLength();
  double endRange   = sStart + getWindingLength();

  if (startRange <= _helixStart) startRange = _helixStart;
  if (endRange >= _helixEnd) endRange = _helixEnd;

  ///###  int nSteps = 100;
  int nSteps = 10;
  double stepWidth = (endRange - startRange)/((double)nSteps);

  double sOfMin = 0;
  double distMinSQ = DBL_MAX;
  double s = startRange;
  LCVector3D x;
  for (int i = 0; i <= nSteps; i++)
    {
      x = getPosition(s);
      double distsq = (x-position).mag2() ;
// ###### vector
      if (distsq < distMinSQ)
	{
	  distMinSQ = distsq;
	  sOfMin = s;
	}
      s += stepWidth;
    }

  double minStepWidth = 0.000001;
  while (stepWidth > minStepWidth)
    {
      stepWidth /= 10;

      x = getPosition(sOfMin + stepWidth);
      double lp = (x-position).mag2();
      x = getPosition(sOfMin - stepWidth);
      double lm = (x-position).mag2();

      double upOrDown = 0;
      double lastDistance = 0;
      nSteps = 10 ;
      if (distMinSQ <= lp && distMinSQ <= lm)
	{
	  nSteps = 0 ;
	  lastDistance = distMinSQ;
	}
      else if (lm<lp)
	{
	  upOrDown = -1;
	  lastDistance = lm;
	}
      else if (lp<lm)
	{
	  upOrDown = 1;
	  lastDistance = lp;
	}
      else
	{
	  //	  cout << "Da laeuft was falsch!!! " 
	  //	       << " distMinSQ " << distMinSQ << " lm " << lm << " lp " << lp << endl;
	}

      double step = 0, l = 0, lastStep = sOfMin;
      for (int i = 1; i<=nSteps;i++)
	{
	  step = sOfMin + (upOrDown*stepWidth*(double)i);
	  x = getPosition(step);
	  l = (x-position).mag2();
	  if (l <= lastDistance) 
	    {
	      lastStep = step;
	      lastDistance = l;
	    }
	}
      sOfMin = lastStep ;
      distMinSQ = lastDistance;
    }
  return sOfMin;
}

double SimpleHelix::getIntersectionWithPlane( LCPlane3D p, 
					      bool& pointExists) const
{

  // Calculation of the region of s (parameter of Helix), where the helix
  // interacts with the plane. This is done by two lines parallel to the 
  // axis of the helix (in this parametrisation the z-axis) in the distance 
  // of the helix radius. Both lines on opposite sides of the helix. 
  // One on the side nearest to the plane, the other exactly on the 
  // other side, away from the plane with respect to the Origen.

  LCVector3D helixCentre( getCentreX(), getCentreY(), 0.);
  LCVector3D normalXY = p.normal();
  normalXY.setZ(0.);
  // If plane is perpendicular to helix axis, set normalXY to a value
  // perpendicular to the helix axis to enable intersection calculation.
  if (normalXY.x() == 0. && normalXY.y() == 0. ) normalXY.set(1.,0.,0.);
  normalXY = ( normalXY.unit() )/_omega;
  LCVector3D lineDirection(0.,0.,1.);

  LCLine3D frontLine( (helixCentre + normalXY) , lineDirection);
  LCLine3D backLine( (helixCentre - normalXY) , lineDirection);
  bool parallelBack, parallelFront;
  double zStart = ( frontLine.position( frontLine.intersectionWithPlane(p,parallelFront) ) ).z();
  double zEnd = ( backLine.position( backLine.intersectionWithPlane(p,parallelBack) ) ).z();
  if ( !(parallelBack && parallelFront) )
    { // plane is parallel to helix axis
      if ( fabs( p.distance(helixCentre) ) > (1/_omega) )
	{ // Helix never hits plane!
	  pointExists = false ;
	  return 0;
	}
      else
	{
	  pointExists = true ;
	  zStart = -DBL_MAX ;
	  zEnd   =  DBL_MAX ;
	}
    }

  double sStart = ( zStart - _reference.z() - _z0)
    * sqrt(1 + _tanLambda*_tanLambda) / _tanLambda;
  double sEnd = ( zEnd - _reference.z() - _z0)
    * sqrt(1 + _tanLambda*_tanLambda) / _tanLambda;

  if (sStart > sEnd)
    { // rong order for Start and Endpoint
      double temp = sEnd;
      sEnd = sStart;
      sStart = temp;
    }

  // produce an artificial gap between sStart and sEnd
  if (sStart > -DBL_MAX && sStart < DBL_MAX ) sStart -= 0.00001;
  if (sEnd > -DBL_MAX && sEnd < DBL_MAX ) sEnd += 0.00001;

  if ( (sStart < 0.) && (sEnd < 0.) )
    { // Intersection is in backwards direction
      pointExists = false ;
      return 0;
    }
  else if ( (sStart < 0.) && (sEnd > 0.) )
    { // intersection region starts in backwards direction
      sStart = 0;
    }

  double s = sStart, sOld = -DBL_MAX; 
  double epsilon = 0.0000001;

  while ( (s-sOld) > epsilon )
    {
      double d = fabs( p.distance( getPosition(s) ) ) ;
      sOld = s;
      s += d;
      if (s > sEnd)
	{ // Problem! no intersection !
	  std::cout << "ERROR: No intersection found!!!" << std::endl;
	  pointExists = false ;
	  return 0;
	}
    }

  pointExists = true ;
  return s ;
}

double SimpleHelix::getIntersectionWithCylinder(const LCCylinder & cylinder,
						bool & pointExists) const  
{
  LCVector3D helixCentre( getCentreX(), getCentreY(), 0.);
  LCVector3D axisDirection(0.,0.,1.);
  LCLine3D axisLine(helixCentre , axisDirection);

  LCVector3D middlePoint =
    (cylinder.startPoint() + cylinder.endPoint())/2.;
  double minDistance = sqrt( 0.25*cylinder.length()*cylinder.length()
                             + cylinder.radius()*cylinder.radius() ); 
                              

  if ( axisLine.distance(middlePoint) > (minDistance+getRadius()) )
    {
      pointExists = false ;
      return 0 ;
    }

  double sProject = axisLine.projectPoint(middlePoint);
  double sStart = sProject - minDistance ;
  double sEnd   = sProject + minDistance ;

  double sHelixStart = ( (axisLine.position(sStart)).z() - _reference.z() - _z0)
    * sqrt(1 + _tanLambda*_tanLambda) / _tanLambda;
  double sHelixEnd = ( (axisLine.position(sEnd)).z() - _reference.z() - _z0)
    * sqrt(1 + _tanLambda*_tanLambda) / _tanLambda;

  if (sHelixStart > sHelixEnd)
    {
      double temp = sHelixStart;
      sHelixStart = sHelixEnd;
      sHelixEnd = temp;
    }

  if ( (sHelixStart < 0.) && (sHelixEnd < 0.) )
    { // Intersection is in backwards direction
      pointExists = false ;
      return 0;
    }
  else if ( (sHelixStart < 0.) && (sHelixEnd > 0.) )
    { // intersection region starts in backwards direction
      sHelixStart = 0;
    }

  double s = sHelixStart, sOld = -DBL_MAX;
  double epsilon = 0.0000001;

  while ( (s-sOld) > epsilon )
    {
      double d = fabs( cylinder.distance( getPosition(s) ) ) ;
      sOld = s;
      if (s > sHelixEnd)
        { // Problem! no intersection !
          pointExists = false ;
          return 0;
        }
      s += d;
    }

  pointExists = true ;
  return s ;
}

double SimpleHelix::getStart() const
{
  return _helixStart;
}

double SimpleHelix::getEnd() const
{
  return _helixEnd;
}

bool SimpleHelix::setStart(double s)
{
  if (s <= _helixEnd)
    {
      _helixStart = s;
      return true;
    }

  return false;
}

bool SimpleHelix::setEnd(double s)
{
  if (s >= _helixStart)
    {
      _helixEnd = s;
      return true;
    }

  return false;
}

bool SimpleHelix::setStartEnd(double start, double end)
{
  if (start <= end)
    {
      _helixStart = start;
      _helixEnd = end;
      return true;
    }

  return false;
}

void SimpleHelix::printProperties()
{
  using namespace std;
  int colWidth = 10;
  cout << "*******************************************************************************" << endl;
  cout << "*" << endl; 
  cout << "* Helix Parameters:" << endl; 
  cout << "* -----------------" << endl; 
  cout << "*" << endl; 
  cout << "* d_0             " << setw(colWidth) << _d0 << endl; 
  cout << "* phi_0           " << setw(colWidth) << _phi0 << endl; 
  cout << "* Omega           " << setw(colWidth) << _omega
       << " (Radius: " << (1/_omega) << ")" << endl; 
  cout << "* z_0             " << setw(colWidth) << _z0 << endl; 
  cout << "* tan lambda      " << setw(colWidth) << _tanLambda << endl; 
  cout << "*" << endl; 
  cout << "* Reference point " 
       << setw(colWidth) << _reference.x() << " "
       << setw(colWidth) << _reference.y() << " "
       << setw(colWidth) << _reference.z() << endl;
  cout << "* pitch           " << setw(colWidth) << getPitch() << endl; 
  cout << "* Winding length  " << setw(colWidth) << getWindingLength() << endl; 
  cout << "* Centre of arc   " 
       << setw(colWidth) << getCentreX() << " "
       << setw(colWidth) << getCentreY() << endl;
  cout << "* Begin and end   " 
       << setw(colWidth) << getStart() 
       << setw(colWidth) << getEnd() << endl; 
  cout << "*" << endl; 
  cout << "*******************************************************************************" << endl;
}

double SimpleHelix::getPitch()
{
  return 2*_pi*_tanLambda/fabs(_omega);
}

double SimpleHelix::getRadius() const 
{
  return 1/_omega;
}

