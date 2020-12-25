#ifndef SecondAlg_h
#define SecondAlg_h

#include <GaudiKernel/Algorithm.h>
#include <Gaudi/Property.h>

#include "Examples/IFirstSvc.h"

// The second algorithm shows how to invoke the services.

class SecondAlg: public Algorithm {
public:
    SecondAlg(const std::string& name, ISvcLocator* pSvcLocator);

    StatusCode initialize() override;
    StatusCode execute() override;
    StatusCode finalize() override;

private:
    SmartIF<IFirstSvc> m_firstsvc;
};


#endif
