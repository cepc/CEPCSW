#ifndef __ILDSITCYLINDERKALDETECTOR__
#define __ILDSITCYLINDERKALDETECTOR__

/** SIT Cylinder based detector to be used for ILD DBD studies when using the old LOI base SIT 
 *
 * @author S.Aplin DESY
 */

#include "kaltest/TVKalDetector.h"

class TNode;

namespace gear{
  class GearMgr ;
}


class ILDSITCylinderKalDetector : public TVKalDetector {
public:
  
  /** Initialize the TPC from GEAR */
  ILDSITCylinderKalDetector( const gear::GearMgr& gearMgr );
  
  
private:
  
  void setupGearGeom( const gear::GearMgr& gearMgr ) ;
  
  unsigned int _nLayers ;
  double _bZ ;
  
  struct SIT_Layer {
    double radius;
    double half_length;
    double senThickness;
    double supThickness;
    
  };
  std::vector<SIT_Layer> _SITgeo;

  
};

#endif
