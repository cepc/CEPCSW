#include "DriftChamberSensDetTool.h"

#include "G4VSensitiveDetector.hh"

#include "DD4hep/Detector.h"

#include "DriftChamberSensitiveDetector.h"
#include "TrackerCombineSensitiveDetector.h"
DECLARE_COMPONENT(DriftChamberSensDetTool)

StatusCode DriftChamberSensDetTool::initialize() {
    StatusCode sc;

    m_geosvc = service<IGeomSvc>("GeomSvc");
    if (!m_geosvc) {
        error() << "Failed to find GeomSvc." << endmsg;
        return StatusCode::FAILURE;
    }

    m_dedx_simtool = ToolHandle<IDedxSimTool>(m_dedx_sim_option.value());
    if (!m_dedx_simtool) {
        error() << "Failed to find dedx simtoo." << endmsg;
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

    if (name == "DriftChamber") {
      if(m_sdTypeOption==0){
	DriftChamberSensitiveDetector* dcsd = new DriftChamberSensitiveDetector(name, *dd4hep_geo);
	dcsd->setDedxSimTool(m_dedx_simtool);
	sd = dcsd;
      }
      else if(m_sdTypeOption==1){
	sd = new TrackerCombineSensitiveDetector(name, *dd4hep_geo);
      }
    }
    else{
      warning() << "The SD " << name << " want to use SD for DriftChamber" << endmsg;
    }

    return sd;
}
