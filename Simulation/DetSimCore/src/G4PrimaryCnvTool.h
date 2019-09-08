#ifndef G4PrimaryCnvTool_h
#define G4PrimaryCnvTool_h

#include "GaudiKernel/AlgTool.h"
#include "DetSimInterface/IG4PrimaryCnvTool.h"

class G4PrimaryCnvTool: public extends<AlgTool, IG4PrimaryCnvTool> {
public:

    using extends::extends;

    bool mutate(G4Event* anEvent) override;

};

#endif
