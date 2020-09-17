#ifndef DriftChamberSensDetTool_h
#define DriftChamberSensDetTool_h

/*
 * DriftChamberSensDetTool is used to create Drift Chamber SD.
 * 
 * It will use DedxSimTool to give the dE/dx value.
 *
 * -- 17 Sept 2020, Tao Lin <lintao@ihep.ac.cn>
 */

#include "GaudiKernel/AlgTool.h"
#include "DetSimInterface/ISensDetTool.h"
#include "DetInterface/IGeoSvc.h"

class DriftChamberSensDetTool: public extends<AlgTool, ISensDetTool> {

public:

    using extends::extends;

    /// Overriding initialize and finalize
    StatusCode initialize() override;
    StatusCode finalize() override;

    /// Override ISensDetTool
    virtual G4VSensitiveDetector* createSD(const std::string& name) override;

private:

    // in order to initialize SD, we need to get the lcdd()
    SmartIF<IGeoSvc> m_geosvc;

};

#endif
