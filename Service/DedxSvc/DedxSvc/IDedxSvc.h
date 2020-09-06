#ifndef I_Dedx_SVC_H
#define I_Dedx_SVC_H

#include "GaudiKernel/IService.h"

#include "G4Step.hh"

class IDedxSvc: virtual public IService {
public:
    DeclareInterfaceID(IDedxSvc, 0, 1); // major/minor version
    
    virtual ~IDedxSvc() = default;
    virtual float pred(const G4Step* aStep)=0 ;

};


#endif
