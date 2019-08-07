#ifndef IFirstSvc_h
#define IFirstSvc_h

#include "GaudiKernel/IService.h"

class IFirstSvc: virtual public IInterface {
public:
    DeclareInterfaceID(IFirstSvc, 0, 1); // major/minor version

    virtual void shoot() = 0;

};

#endif
