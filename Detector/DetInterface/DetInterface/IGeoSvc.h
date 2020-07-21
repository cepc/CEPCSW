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
#include "DDRec/DetectorData.h"
#include <map>

namespace dd4hep {
  class Detector;
  class DetElement;
  // if not include DDRec/DetectorData.h, can not compile, why? 
  //namespace rec{
    //struct StructExtension;
    //struct ZPlanarStruct;
    //typedef StructExtension<ZPlanarStruct> ZPlanarData;
    //struct ConicalSupportStruct;
    //typedef StructExtension<ConicalSupportStruct> ConicalSupportData;
}

class StructExtension;

namespace gear{
  class ZPlanarParametersImpl;
  class GearParametersImpl;
}
class TMaterial;
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

  // obsolete parameter format, will remove once StructExtension<> validated
  virtual const gear::ZPlanarParametersImpl*  getVXDParameters() = 0;

  virtual const dd4hep::rec::ZPlanarData* getVXDData() = 0;
  virtual const dd4hep::rec::ConicalSupportData* getBeamPipeData() =0;

  virtual const std::map<std::string,double>& getDetParameters(std::string s) = 0;
  virtual const double getDetParameter(std::string set_name, std::string par_name) = 0;
  virtual TMaterial* getMaterial(std::string s) = 0;
  virtual ~IGeoSvc() {}
};

#endif  // IGEOSVC_H
