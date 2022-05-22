#include "GenericTrackerSensDetTool.h"

#include "G4VSensitiveDetector.hh"

#include "DD4hep/Detector.h"

#include "GenericTrackerSensitiveDetector.h"

#include "CLHEP/Units/SystemOfUnits.h"

DECLARE_COMPONENT(GenericTrackerSensDetTool)

StatusCode GenericTrackerSensDetTool::initialize() {
  StatusCode sc;
  debug() << "initialize() " << endmsg;
  
  m_geosvc = service<IGeomSvc>("GeomSvc");
  if (!m_geosvc) {
    error() << "Failed to find GeomSvc." << endmsg;
    return StatusCode::FAILURE;
  }
  
  return AlgTool::initialize();
}

StatusCode GenericTrackerSensDetTool::finalize() {
  StatusCode sc;
  
  return sc;
}

G4VSensitiveDetector* GenericTrackerSensDetTool::createSD(const std::string& name) {
  debug() << "createSD for " << name << endmsg;
  
  dd4hep::Detector* dd4hep_geo = m_geosvc->lcdd();

  G4VSensitiveDetector* sd = new GenericTrackerSensitiveDetector(name, *dd4hep_geo);

  return sd;
}
