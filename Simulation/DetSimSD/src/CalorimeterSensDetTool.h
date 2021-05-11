#ifndef CalorimeterSensDetTool_h
#define CalorimeterSensDetTool_h

/*
 * CalorimeterSensDetTool is used to create the Calorimeter SD.
 *
 * -- 12 June 2020, Tao Lin <lintao@ihep.ac.cn>
 */

#include "GaudiKernel/AlgTool.h"
#include "DetSimInterface/ISensDetTool.h"
#include "DetInterface/IGeomSvc.h"

class CalorimeterSensDetTool: public extends<AlgTool, ISensDetTool> {

public:

    using extends::extends;

    /// Overriding initialize and finalize
    StatusCode initialize() override;
    StatusCode finalize() override;

    /// Override ISensDetTool
    virtual G4VSensitiveDetector* createSD(const std::string& name) override;

private:

    // in order to initialize SD, we need to get the lcdd()
    SmartIF<IGeomSvc> m_geosvc;

};

#endif
