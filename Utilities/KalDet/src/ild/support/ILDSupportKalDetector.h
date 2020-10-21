#ifndef __ILDSUPPORTDETECTOR__
#define __ILDSUPPORTDETECTOR__

/** Support Material to be used for ILD DBD studies 
 *
 * @author S.Aplin DESY
 */

#include "kaltest/TVKalDetector.h"

class TNode;

namespace gear{
  class GearMgr ;
}
class IGeomSvc;

class ILDCylinderMeasLayer;

class ILDSupportKalDetector : public TVKalDetector {
public:
  
  /** Initialize the support structures from GEAR */
  ILDSupportKalDetector( const gear::GearMgr& gearMgr, IGeomSvc* geoSvc );
  
  /** Returns the special layer inside the Beam Pipe used for propagation to the IP */
  ILDCylinderMeasLayer* getIPLayer() { return _ipLayer; }
  
private:
  
  ILDCylinderMeasLayer* _ipLayer;
  
  std::vector<std::string> _surface_names;
  
};

#endif
