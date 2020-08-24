#ifndef MiniVectorHit01_h
#define MiniVectorHit01_h

#include "ILDImpl/IMiniVector.h"

namespace KiTrackMarlin{
  /** A class for mini-vectors in the VXD - SIT system (the 01 is just for historical reasons and may be renamed)
   */   
  class MiniVectorHit01 : public IMiniVector{
  public:
    
    MiniVectorHit01( MiniVector* miniVector , const SectorSystemVXD* const sectorSystemVXD );
  };
}
#endif

