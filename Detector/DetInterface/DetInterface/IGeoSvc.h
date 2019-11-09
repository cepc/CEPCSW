//
//  IGeoSvc.h
//
//  Based on FCCSW with some modification.
//  In the design, the geometry shoud only depends on DD4hep. 
//  
//  -- Tao Lin, 2019/11/09
//

#ifndef IGEOSVC_H
#define IGEOSVC_H

#include "GaudiKernel/IService.h"

namespace dd4hep {
class Detector;
class DetElement;
}

// class G4VUserDetectorConstruction;

class GAUDI_API IGeoSvc : virtual public IService {

public:
  /// InterfaceID
  DeclareInterfaceID(IGeoSvc, 1, 0);
  // receive DD4hep Geometry
  virtual dd4hep::DetElement getDD4HepGeo() = 0;
  virtual dd4hep::Detector* lcdd() = 0;
  // receive Geant4 Geometry
  // virtual G4VUserDetectorConstruction* getGeant4Geo() = 0;

  virtual ~IGeoSvc() {}
};

#endif  // IGEOSVC_H
