#ifndef I_GEAR_SVC_H
#define I_GEAR_SVC_H

#include "GaudiKernel/IService.h"
#include "gear/GearMgr.h"

// IGearSvc is the interface between Gaudi and GEAR.

class IGearSvc: virtual public IService {
public:
    DeclareInterfaceID(IGearSvc, 0, 1); // major/minor version
    
    virtual ~IGearSvc() = default;

    // Get the GEAR Manager
    virtual gear::GearMgr* getGearMgr() = 0;

};


#endif
