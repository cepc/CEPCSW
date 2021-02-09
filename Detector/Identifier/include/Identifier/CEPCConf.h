/* CEPC conference description ID */
#ifndef CEPCConf_H
#define CEPCConf_H

#include <string>

namespace CEPCConf{
  struct DetID{ // compatible for old codes from Marlin, will convert from compact file, default initial values here
    static const int NOTUSED = 0;
    static const int VXD     = 1;
    static const int SIT     = 2;
    static const int FTD     = 3;
    static const int TPC     = 4;
    static const int SET     = 5;
    static const int ETD     = 6;
    static const int DC      = 7; 
    
    static const int ECAL        = 20;
    static const int ECAL_PLUG   = 21;
    static const int HCAL        = 22;
    static const int HCAL_RING   = 23;
    static const int LCAL        = 24;
    static const int BCAL        = 25;
    static const int LHCAL       = 26;
    static const int YOKE        = 27;
    static const int COIL        = 28;
    static const int ECAL_ENDCAP = 29;
    static const int HCAL_ENDCAP = 30;
    static const int YOKE_ENDCAP = 31;
    
    static const int bwd    = -1;
    static const int barrel =  0;
    static const int fwd    = +1;
  };
  
  struct TrkHitTypeBit{
    static const int ONE_DIMENSIONAL      = 29;
    static const int COMPOSITE_SPACEPOINT = 30;
    static const int PLANAR               = 3; // 3 is compatible with old tracking codes, 31 or 28 is better in future to modify uniformly
  };
  
  struct TrkHitQualityBit{
    static const int USED_IN_FIT          = 30;
    static const int USED_IN_TRACK        = 29;
    static const int DOUBLE_HIT_CANDIDATE = 28;
    static const int GOOD                 = 27;
  };
}
#endif
