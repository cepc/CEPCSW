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
    CaloSensitiveDetector* sd = new CaloSensitiveDetector(name, *dd4hep_geo, is_merge_enabled);
    debug() << name << " set to merge true/false = " << is_merge_enabled << endmsg;

    auto sens = dd4hep_geo->sensitiveDetector(name);
    std::string typ = sens.type();
    if(typ=="scintillator"&&m_applyBirksLaw) sd->ApplyBirksLaw();

    return sd;
}


