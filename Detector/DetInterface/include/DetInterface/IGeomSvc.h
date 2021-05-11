//
//  IGeomSvc.h
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
    namespace DDSegmentation {
        class BitFieldCoder;
    }
  class Detector;
  class DetElement;

}

class StructExtension;

namespace gear{
  class ZPlanarParametersImpl;
  class GearParametersImpl;
}
class TMaterial;
// class G4VUserDetectorConstruction;

class GAUDI_API IGeomSvc : virtual public IService {
public:
  typedef dd4hep::DDSegmentation::BitFieldCoder Decoder;
public:
  /// InterfaceID
  DeclareInterfaceID(IGeomSvc, 1, 0);
  // receive DD4hep Geometry
  virtual dd4hep::DetElement getDD4HepGeo() = 0;
  virtual dd4hep::Detector* lcdd() = 0;
  // receive Geant4 Geometry
  // virtual G4VUserDetectorConstruction* getGeant4Geo() = 0;

  // short cut to retrieve the Decoder according to the Readout name
  virtual Decoder* getDecoder(const std::string& readout_name) = 0;

  // obsolete parameter format, will remove once StructExtension<> validated
  virtual const gear::ZPlanarParametersImpl*  getVXDParameters() = 0;

  virtual const dd4hep::rec::ZPlanarData* getVXDData() = 0;
  virtual const dd4hep::rec::ConicalSupportData* getBeamPipeData() =0;

  virtual const std::map<std::string,double>& getDetParameters(std::string s) = 0;
  virtual double getDetParameter(std::string set_name, std::string par_name) = 0;
  virtual TMaterial* getMaterial(std::string s) = 0;

  virtual ~IGeomSvc() {}
};

#endif  // IGEOSVC_H
