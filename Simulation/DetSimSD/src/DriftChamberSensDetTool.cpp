#include "DriftChamberSensDetTool.h"

#include "G4VSensitiveDetector.hh"

#include "DD4hep/Detector.h"

DECLARE_COMPONENT(DriftChamberSensDetTool);

StatusCode DriftChamberSensDetTool::initialize() {
    StatusCode sc;

    m_geosvc = service<IGeoSvc>("GeoSvc");
    if (!m_geosvc) {
        error() << "Failed to find GeoSvc." << endmsg;
        return StatusCode::FAILURE;
    }

    return sc;
}

StatusCode DriftChamberSensDetTool::finalize() {
    StatusCode sc;
    
    return sc;
}

G4VSensitiveDetector*
DriftChamberSensDetTool::createSD(const std::string& name) {
    dd4hep::Detector* dd4hep_geo = m_geosvc->lcdd();

    G4VSensitiveDetector* sd = nullptr;

    return sd;
}
