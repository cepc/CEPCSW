#ifndef VXDHit01_h
#define VXDHit01_h

#include "ILDImpl/IVXDHit.h"

namespace KiTrackMarlin{
  /** A class for hits in the VXD (the 01 is just for historical reasons and may be renamed)
   * 
   * - The side is according to CellID0.
   * - Layer is set according to CellID0 +1 (so we can use layer 0 for the IP)
   * - Module is set according to CellID0.
   * - Sensor is set according to CellID0 -1. (because currently sensors of the VXD start with 1 in the CellID0, if this changes, this has to be modified)
   */   
  class VXDHit01 : public IVXDHit{
  public:
    
    VXDHit01( edm4hep::TrackerHit* trackerHit , const SectorSystemVXD* const sectorSystemVXD );
  };
  //void setSectorisationInPhi(int PhiSectors);
  //void setSectorisationInTheta(int ThetaSectors);
}
#endif

