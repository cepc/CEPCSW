#include "CalorimeterSensDetTool.h"

#include "G4VSensitiveDetector.hh"

#include "DetSimSD/CaloSensitiveDetector.h"

#include "DD4hep/Detector.h"

DECLARE_COMPONENT(CalorimeterSensDetTool);

StatusCode
CalorimeterSensDetTool::initialize() {
    StatusCode sc;


    m_geosvc = service<IGeomSvc>("GeomSvc");
    if (!m_geosvc) {
        error() << "Failed to find GeomSvc." << endmsg;
        return StatusCode::FAILURE;
    }


    return sc;
}

StatusCode
CalorimeterSensDetTool::finalize() {
    StatusCode sc;

    return sc;
}

G4VSensitiveDetector*
CalorimeterSensDetTool::createSD(const std::string& name) {

    dd4hep::Detector* dd4hep_geo = m_geosvc->lcdd();

    bool is_merge_enabled = true;
    for(auto cal_name : m_listCalsMergeDisable){
      if(cal_name==name){
	is_merge_enabled = false;
	break;
      }
    }
    G4VSensitiveDetector* sd = new CaloSensitiveDetector(name, *dd4hep_geo, is_merge_enabled);
    warning() << name << " set to merge true/false = " << is_merge_enabled << endmsg;

    return sd;
}


