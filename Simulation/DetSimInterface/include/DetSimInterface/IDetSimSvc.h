#ifndef IDetSimSvc_h
#define IDetSimSvc_h

#include "GaudiKernel/IService.h"

// IDetSimSvc is the interface between Gaudi and Geant4.
// All the initialization of Run Manager (RM) should be done here, including:
//   * Detector Construction
//   * Physics List
//   * Primary Generator Action
//   * User Actions
// Then, the real simulation should be also done by this service via Run Manager.
//
// Note, to decouple the Gaudi and Geant4, we keep all these classes still derived from
// the original Geant4's base classes, while using Gaudi tools to manage these objects.

class G4RunManager;

class IDetSimSvc: virtual public IInterface {
public:
    DeclareInterfaceID(IDetSimSvc, 0, 1); // major/minor version
    
    virtual ~IDetSimSvc() = 0;

    // Get the Run Manager
    virtual G4RunManager* getRM() = 0;

    // Control the run manager directly.
    virtual StatusCode initializeRM() = 0;
    virtual StatusCode simulateEvent(int i_event) = 0;
    virtual StatusCode finalizeRM() = 0;

};


#endif
