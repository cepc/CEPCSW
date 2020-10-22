#ifndef LCLine3D_H
#define LCLine3D_H 1

// #include "CLHEP/Vector/ThreeVector.h"
#include <DataHelper/LCGeometryTypes.h>
#include <DataHelper/LCPlane3D.h>

/** Definition of a LCLine3D describing a geometrical line in 3D space.
 *  @author T.Kraemer, DESY
 *  @version $Id: LCLine3D.h,v 1.8 2006-11-03 16:43:13 tkraemer Exp $
 */

class LCLine3D {

public:

  /** Standard constructor:
   * Initializes a line along the x-axis.
   */
  LCLine3D();

  /**
   * Constructor from a point and a direction.
   * @param point Point is a point of the line
   * @param direction Direction is the directional vector of the line.
   */
  LCLine3D(const LCVector3D & point, const LCVector3D & direction) ;

  /**
   * Constructor from a point and a direction.
   * @param point Point is a point of the line
   * @param direction Direction is the directional vector of the line.
   * @param reference reference point of the line. 
   */
  LCLine3D(const LCVector3D & point, 
	   const LCVector3D & direction, 
	   const LCVector3D & reference) ;

  /**
   * Constructor using the canonical parameterization. 
   * @param d0 d0 is the point of closest approach in the xy plane.
   * @param phi0 phi0 is the angle in  the xy plane.
   * @param z0 z0 is the z coordinate of the point of closest approach.
   * @param tanLambda tanLambda is the angle of with respect to the xy plane.
   */
  LCLine3D(double d0, double phi0, double z0, double tanLambda) ;

  /**
   * Constructor using the canonical parameterization. 
   * @param d0 d0 is the point of closest approach in the xy plane.
   * @param phi0 phi0 is the angle in  the xy plane.
   * @param z0 z0 is the z coordinate of the point of closest approach.
   * @param tanLambda tanLambda is the angle of with respect to the xy plane.
   * @param reference reference point of the line. 
   */
  LCLine3D(double d0, double phi0, double z0, double tanLambda,
	   const LCVector3D & reference) ;

  /** Copy constructor.
   * @param line line is an other LCLine3D.
   */
  LCLine3D(const LCLine3D & line) ;

  /**
   * Destructor. */
  ~LCLine3D() {}

  /**
   * set the Parameters for a line using a point and a direction.
   * @param point Point is a point of the line
   * @param direction Direction is the directional vector of the line.
   * @param reference reference point of the line. 
   */
  bool set(const LCVector3D & point, 
	   const LCVector3D & direction, 
	   const LCVector3D & reference) ;

  /**
   * Set the Parameters of a line using the canonical parameterization. 
   * @param d0 d0 is the point of closest approach in the xy plane.
   * @param phi0 phi0 is the angle in  the xy plane.
   * @param z0 z0 is the z coordinate of the point of closest approach.
   * @param tanLambda tanLambda is the angle of with respect to the xy plane.
   * @param reference reference point of the line. 
   */
  bool set(double d0, double phi0, double z0, double tanLambda,
	   const LCVector3D & reference) ;

  /**
   * Assignment. */
  LCLine3D & operator=(const LCLine3D & rhs) ;

  /**
   * Position is the point of the line after a distance s. 
   * Is is given with respect to the point of closes approach to the origen of
   * the coordinate system. 
   * @param s s is the path length along the line */
  LCVector3D position(const double s = 0) const ;

  /** Direction of the line 
   */
  LCVector3D direction() const ;

  /**
   * Distance of a point to the line. 
   * @param point point is a point in space
   */
  double distance(const LCVector3D & point) const ;  

  /**
   * Projection of a point on to the line. 
   * @param point point is a point in space.
   */
  double projectPoint(const LCVector3D & point) const ;

  /**
   * Test for equality. */
  bool operator==(const LCLine3D & rhs) const ; 

  /**
   * Test for inequality. */
  bool operator!=(const LCLine3D & rhs) const ; 

  /** Pathlength at closest intersection point with plane - undefined
   *  if pointExists==false.
   */
  double intersectionWithPlane(const LCPlane3D plane, bool& pointExists) const ;

protected:

  LCVector3D _point;
  LCVector3D _direction;
  LCVector3D _reference;
};

std::ostream & operator << (std::ostream &os, const LCLine3D &l) ;

#endif /* ifndef LCLine3D_H */
