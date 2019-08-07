#ifndef FirstSvc_h
#define FirstSvc_h

#include "Examples/IFirstSvc.h"
#include <GaudiKernel/Service.h>

class FirstSvc: public extends<Service, IFirstSvc> {
public:
    using extends::extends;

    StatusCode initialize() override;
    StatusCode finalize() override;

    void shoot() override;
};

#endif
