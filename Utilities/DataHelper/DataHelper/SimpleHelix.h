#ifndef SimpleHelix_H
#define SimpleHelix_H 1

#include "DataHelper/Trajectory.h"

/** Simple helix trajectory.
 *  @author T.Kraemer, DESY
 *  @version $Id: SimpleHelix.h,v 1.7 2007-06-20 18:47:25 samson Exp $
 */

class SimpleHelix : public Trajectory {

public:

  virtual ~SimpleHelix() {} 
 
  /** Construct Helix from canonical parameters.
   */
  SimpleHelix( double d0, double phi0, double omega,
	       double z0, double tanLambda, 
	       LCVector3D referencePoint, LCErrorMatrix* errors=0) ;
  
  /** Position at path length s - s==0 corresponds to P.C.A to the origin.
   *  @param s      path length
   *  @param errors return argument - not computed if NULL
   */
  virtual LCVector3D getPosition(double s, LCErrorMatrix* errors=0) const ;
  
  /** Direction at path length s, i.e. (dx/ds,dy/ds,dz/ds) 
   *  @param s      path length
   *  @param errors return argument - not computed if NULL
   */
  virtual LCVector3D getDirection(double s,  LCErrorMatrix* errors=0) const ;
  
  /** Full covariance Matrix of x,y,z,px,py,pz   
   *  @param s      path length
   */
  virtual LCErrorMatrix getCovarianceMatrix( double s) const ;

  /** Pathlength at point on trajectory closest to given position.  
   *  In order to get the distance use for example:  <br>  
   *     LCVector3D pt = t.getPosition( t.getPathAtClosestPoint( p ) ) ; <br>
   *     double d = LCVector3D( pt - p ).mag()  ; <br> 
   */
  virtual double getPathAt(const LCVector3D position ) const ;
  
  /** get the radius of the helix. 
   */
  virtual double getRadius() const ;

  /*----------------------------------------------------------------------*/
  
  /** Pathlength at closest intersection point with plane - undefined 
   *  if pointExists==false. 
   */
  virtual double getIntersectionWithPlane( LCPlane3D p, bool& pointExists) const  ;
  
  /** Pathlength at closest intersection point with cylinder - undefined 
   *  if pointExists==false. 
   * @param cylinder cylinder object to intersect with
   */
  virtual  double getIntersectionWithCylinder(const LCCylinder & cylinder,
                                              bool & pointExists) const ;

  /** Pathlength at the start and end point of the trajectory. 
   */
  virtual double getStart() const ;
  virtual double getEnd() const ;

  /** Set pathlength at the start and end point of the trajectory. 
   * @param start pathlength of start point 
   * @param end pathlength of end point 
   */
  virtual bool setStart(double s);
  virtual bool setEnd(double s);
  virtual bool setStartEnd(double start, double end);

  virtual void printProperties();

protected:

  SimpleHelix() {} 

  virtual void init() ;

  virtual double getCentreX() const ;
  virtual double getCentreY() const ;
  virtual double getWindingLength() const ;
  virtual double getPitch();
  double _d0;
  double _phi0;
  double _omega;
  double _z0;
  double _tanLambda;

  //  double _Bz;

  double _helixStart;
  double _helixEnd;

  static const char* _names[];

  static const double _a; // = 2.99792458E-4;
  static const double _pi; // = 3.14159265358979323846;

  LCVector3D  _reference; 
  LCErrorMatrix* _errors;
}; // class 



// /** Physical trajectory describing a (charged) particle's  path in a B 
//  *  field and material. 
//  *  @author F.Gaede, DESY
//  *  @version $Id: SimpleHelix.h,v 1.7 2007-06-20 18:47:25 samson Exp $
//  */

// class PhysicalSimpleLine : public SimpleLine{
  
//   /** Particle's momentum at path length s. Implementations will have to  have knowledge
//    *  about the particle type, B-field and material.
//    */
//   virtual LCLorentzVector get4Momentum( double s ) const ;
// }


#endif /* ifndef SimpleLine_H */
