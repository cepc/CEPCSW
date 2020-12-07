#ifndef DummyFastSimG4Tool_h
#define DummyFastSimG4Tool_h

#include "GaudiKernel/AlgTool.h"
#include "DetSimInterface/IFastSimG4Tool.h"

class DummyFastSimG4Tool: public extends<AlgTool, IFastSimG4Tool> {
public:
    using extends::extends;

    StatusCode initialize() override;
    StatusCode finalize() override;

    bool CreateFastSimulationModel() override;

private:

};

#endif
