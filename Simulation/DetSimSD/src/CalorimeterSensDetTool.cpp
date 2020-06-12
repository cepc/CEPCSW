#include "CalorimeterSensDetTool.h"

#include "G4VSensitiveDetector.hh"

DECLARE_COMPONENT(CalorimeterSensDetTool);

StatusCode
CalorimeterSensDetTool::initialize() {
    StatusCode sc;

    return sc;
}

StatusCode
CalorimeterSensDetTool::finalize() {
    StatusCode sc;

    return sc;
}

G4VSensitiveDetector*
CalorimeterSensDetTool::createSD(const std::string& name) {

    return nullptr;
}


