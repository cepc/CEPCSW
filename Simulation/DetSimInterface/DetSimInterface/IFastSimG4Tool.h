#ifndef IFastSimG4Tool_h
#define IFastSimG4Tool_h

// IFastSimG4Tool is to associate the G4Region and Fast simulation model in G4.
// It is recommended to create one fast simulation model in one tool.
// -- Tao Lin <lintao@ihep.ac.cn>, 7 Dec 2020

#include "GaudiKernel/IAlgTool.h"

class IFastSimG4Tool: virtual public IAlgTool {
public:
    DeclareInterfaceID(IFastSimG4Tool, 0, 1);

    virtual ~IFastSimG4Tool() {}

    // Build the association between G4Region and G4 Fast simulation model
    virtual bool CreateFastSimulationModel() = 0;
};

#endif
