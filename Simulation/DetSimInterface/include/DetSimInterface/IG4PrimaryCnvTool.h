#ifndef IG4PrimaryCnvTool_h
#define IG4PrimaryCnvTool_h

// IG4PrimaryCnvTool:
//     convert an event in other formats to an event in G4 format.
// The G4Event object is managed by Geant4.

#include "GaudiKernel/AlgTool.h"

class G4Event;

class IG4PrimaryCnvTool: virtual public IAlgTool {
public:

    DeclareInterfaceID(IG4PrimaryCnvTool, 0, 1);

    virtual ~IG4PrimaryCnvTool() {};

    virtual bool mutate(G4Event* anEvent) = 0;
};

#endif
