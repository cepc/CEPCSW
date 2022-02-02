#ifndef FTDHit01_h
#define FTDHit01_h

#include "ILDImpl/IFTDHit.h"

namespace KiTrackMarlin{
  /** A class for hits in the FTD (the 01 is just for historical reasons and may be renamed)
   * 
   * - The side is according to CellID0.
   * - Layer is set according to CellID0 +1 (so we can use layer 0 for the IP)
   * - Module is set according to CellID0.
   * - Sensor is set according to CellID0 -1. (because currently sensors of the FTD start with 1 in the CellID0, if this changes, this has to be modified)
   */   
  class FTDHit01 : public IFTDHit{
  public:
      
    FTDHit01( edm4hep::TrackerHit trackerHit , const SectorSystemFTD* const sectorSystemFTD );
  };
}
#endif

