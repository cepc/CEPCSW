#ifndef LCCylinder_H
#define LCCylinder_H 1

// #include "CLHEP/Vector/ThreeVector.h"
#include <DataHelper/LCGeometryTypes.h>

/** Definition of a LCCylinder describing a geometrical cylinder in 3D space.
 *  @author T.Kraemer, DESY
 *  @version $Id: LCCylinder.h,v 1.1 2006-10-16 15:23:32 tkraemer Exp $
 */

class LCCylinder {

public:

  /**
   * Constructor from two points and a radius. 
   * @param point1 point1 is the start point of the cylinder axis
   * @param point2 point2 is the end point of the cylinder axis
   * @param radius radius is the radius of the xylinder
   * @param endPlane endPlane switches if cylinder is open or not
   */
  LCCylinder(LCVector3D point1, 
	     LCVector3D point2, 
	     double radius, 
	     bool endPlane = false) ;

  /**
   * Constructor from one point, th Axis of the cylinder and the radius
   * @param radius radius of cylinder.
   * @param point point in the middle of the axis of the xylinder
   * @param axis axis is the orientation of the xylinder axis. 
   *             the length of axis is the half length of the cylinder
   * @param endPlane endPlane switches if cylinder is open or not
   */
  LCCylinder(double radius, LCVector3D point, LCVector3D axis, bool endPlane) ;

  /** Copy constructor.
   * @param cylinder cylinder is an other LCCylinder.
   */
  LCCylinder(const LCCylinder & cylinder) ;

  /**
   * Destructor. */
  ~LCCylinder() {}

  /**
   * Assignment. */
  LCCylinder & operator=(const LCCylinder & rhs) ;

  /**
   * startpoint of cylinder axis */
  LCVector3D startPoint() const ;

  /**
   * end point of cylinder axis */
  LCVector3D endPoint() const ;

  /**
   * orientation of cylinder axis. the return vector is normalised */
  LCVector3D axisDirection() const ;

  /**
   * length of cylinder axis */
  double length() const ;

  /**
   * Radius of cylinder.  */
  double radius() const ;

  /**
   * Distance of a point to the cylinder. 
   * @param point point is a point in space
   */
  double distance(const LCVector3D & point) const ;  

  /**
   * Projection of a point on to the surface of the cylinder. 
   * @param point point is a point in space.
   * @param code code gives an integer code for the type of area the 
   * point gets projected to:
   * 0 : No projection possible (point is not normal above the surface)
   * 1 : prjection hits plane at start point ( startPoint() ) 
   * 2 : prjection hits plane at end point ( endPoint() ) 
   * 3 : projection hits the finite tube of the cylinder
   */
  LCVector3D projectPoint(const LCVector3D & point, int & code) const ;

  /**
   * Checks if a given point is inside the clyinder.
   * @param point point is a point in space.
   */
  bool isInside(const LCVector3D & point) const ;

  /**
   * Test for equality. */
  bool operator==(const LCCylinder & rhs) const ; 

  /**
   * Test for inequality. */
  bool operator!=(const LCCylinder & rhs) const ; 

protected:

  double _radius;
  bool _endPlane;

  LCVector3D _axisSstartPoint;
  LCVector3D _axisEndPoint;

};

#endif /* ifndef LCCylinder_H */
