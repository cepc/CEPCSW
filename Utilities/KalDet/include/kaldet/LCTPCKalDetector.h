#ifndef LCTPCKALDETECTOR_H
#define LCTPCKALDETECTOR_H

#include "kaltest/TVKalDetector.h"

#include "ILDVMeasLayer.h"

namespace gear{
  class GearMgr ;
}

namespace kaldet{

  /**
   * The LCTPC implementation for a TPC which is completely instantiated from GEAR.
   * 
   */
class LCTPCKalDetector : public TVKalDetector {

public:

    LCTPCKalDetector() {};

    /** 
     * The constructor. All information to initialise the TPC is taken from GEAR.
     *
     * The class has been copied from GearTPCKalDetector class and adopted for the use of MarlinTrk
     * You can find comments and necessary information in the original class
     * 
     */
    LCTPCKalDetector(const gear::GearMgr& gearMgr);

    /// The destructor.
    virtual ~LCTPCKalDetector();

};

}// namespace kaldet
#endif //LCTPCKALDETECTOR_H
