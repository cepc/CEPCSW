#ifndef __ILDSETKALDETECTOR__
#define __ILDSETKALDETECTOR__

/** Ladder based SET to be used for ILD DBD studies 
 *
 * @author S.Aplin DESY
 */

#include "kaltest/TVKalDetector.h"

#include "TMath.h"

class TNode;

namespace gear{
  class GearMgr ;
}

class IGeomSvc;

class ILDSETKalDetector : public TVKalDetector {
  
public:
  
  /** Initialize the SET from GEAR */
  ILDSETKalDetector( const gear::GearMgr& gearMgr, IGeomSvc* geoSvc=0 );
  
  
private:
  
  void setupGearGeom( const gear::GearMgr& gearMgr ) ;
  void setupGearGeom( IGeomSvc* geoSvc );
  
  int _nLayers ;
  double _bZ ;

  bool _isStripDetector;
    
  struct SET_Layer {
    int nLadders;
    int nSensorsPerLadder;
    double phi0;
    double dphi;
    double senRMin;
    double supRMin;
    double length;
    double width;
    double offset;
    double senThickness;
    double supThickness;
    double sensorLength;
    double stripAngle;
  };
  std::vector<SET_Layer> _SETgeo;
  
};



#endif
