#ifndef LCGeometryTypes_H
#define LCGeometryTypes_H 1

// #include "CLHEP/Geometry/Point3D.h"
// #include "CLHEP/Geometry/Vector3D.h"
#include "CLHEP/Geometry/Plane3D.h"
#include "CLHEP/Matrix/SymMatrix.h"

#include "CLHEP/Vector/ThreeVector.h"


/** @file Definition of geometry types used in ILC software - currently use CHEP.
 *  @author gaede
 *  @version $Id: LCGeometryTypes.h,v 1.3 2006-10-11 16:03:24 tkraemer Exp $
 */

//using namespace CLHEP ;


// typedef HepGeom::Point3D<double> LCPoint3D ;
// typedef HepGeom::Vector3D<double> LCVector3D ;

// typedef CLHEP::Hep3Vector LCPoint3D ;
typedef CLHEP::Hep3Vector LCVector3D ;


typedef CLHEP::HepSymMatrix LCErrorMatrix ;

// typedef HepGeom::Plane3D<double> LCPlane3D ;

typedef CLHEP::HepLorentzVector LCLorentzVector ;






#endif /* ifndef LCGeometryTypes_H */
