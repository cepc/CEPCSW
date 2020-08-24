#ifndef __ILDFTDDISCBASEDDETECTOR__
#define __ILDFTDDISCBASEDDETECTOR__

/** Disk based version of the FTD alla LOI
*
* @author S.Aplin DESY
*/

#include "kaltest/TVKalDetector.h"

class TNode;

namespace gear{
  class GearMgr ;
}


class ILDFTDDiscBasedKalDetector : public TVKalDetector {
public:
  
  /** Initialize the FTD from GEAR */
  ILDFTDDiscBasedKalDetector( const gear::GearMgr& gearMgr );
  
  
private:
  
  void setupGearGeom( const gear::GearMgr& gearMgr ) ;
  
  int _nDisks ;
  double _bZ ;
  
  struct FTD_Disk {
    double rInner;
    double rOuter;
    double senThickness;
    double supThickness;
    double zPos;
    
  };
  std::vector<FTD_Disk> _FTDgeo;
  
};

#endif
