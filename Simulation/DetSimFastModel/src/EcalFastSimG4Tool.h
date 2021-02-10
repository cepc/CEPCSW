#ifndef EcalFastSimG4Tool_h
#define EcalFastSimG4Tool_h

#include "GaudiKernel/AlgTool.h"
#include "DetSimInterface/IFastSimG4Tool.h"

class EcalFastSimG4Tool: public extends<AlgTool, IFastSimG4Tool> {
public:
    using extends::extends;

    StatusCode initialize() override;
    StatusCode finalize() override;

    bool CreateFastSimulationModel() override;

private:
    // the regions will be associated with the fast sim model
    Gaudi::Property<std::vector<std::string>> m_regions{this, "Regions"};
};

#endif
