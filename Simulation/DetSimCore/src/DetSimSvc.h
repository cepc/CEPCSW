#ifndef DetSimSvc_h
#define DetSimSvc_h

#include "DetSimInterface/IDetSimSvc.h"
#include <GaudiKernel/Service.h>

class DetSimSvc: public extends<Service, IDetSimSvc> {
public:

    DetSimSvc(const std::string& name, ISvcLocator* svc );
    ~DetSimSvc();

    // Get the Run Manager
    G4RunManager* getRM() override;

    // Control the run manager directly.
    StatusCode initializeRM() override;
    StatusCode simulateEvent(int i_event) override;
    StatusCode finalizeRM() override;

    StatusCode initialize() override;
    StatusCode finalize() override;

private:
    G4RunManager* m_runmgr;

};


#endif
