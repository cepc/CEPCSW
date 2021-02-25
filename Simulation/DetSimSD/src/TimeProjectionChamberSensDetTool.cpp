#include "TimeProjectionChamberSensDetTool.h"

#include "G4VSensitiveDetector.hh"

#include "DD4hep/Detector.h"

#include "TimeProjectionChamberSensitiveDetector.h"

#include "CLHEP/Units/SystemOfUnits.h"

DECLARE_COMPONENT(TimeProjectionChamberSensDetTool)

StatusCode TimeProjectionChamberSensDetTool::initialize() {
  StatusCode sc;
  debug() << "initialize() " << endmsg;
  debug() << "TypeOption = " << m_sdTypeOption.value() << endmsg;
  
  m_geosvc = service<IGeomSvc>("GeomSvc");
  if (!m_geosvc) {
    error() << "Failed to find GeomSvc." << endmsg;
    return StatusCode::FAILURE;
  }
  
  return AlgTool::initialize();
}

StatusCode TimeProjectionChamberSensDetTool::finalize() {
  StatusCode sc;
  
  return sc;
}

G4VSensitiveDetector* TimeProjectionChamberSensDetTool::createSD(const std::string& name) {
  debug() << "createSD for TPC" << endmsg;
  
  dd4hep::Detector* dd4hep_geo = m_geosvc->lcdd();

  G4VSensitiveDetector* sd = nullptr;
  if (name == "TPC") {
    if(m_sdTypeOption==1/*FIXME: const defined by CEPC whole description*/){
      TimeProjectionChamberSensitiveDetector* tpcsd = new TimeProjectionChamberSensitiveDetector(name, *dd4hep_geo);
      // switch units to Geant4's, Geant4 uses CLHEP units in fact
      tpcsd->setThreshold(m_threshold/dd4hep::eV*CLHEP::eV); //FIXME: let SD get threshold from dd4hep_geo object
      tpcsd->setSameStepLimit(m_sameStepLimit);
      tpcsd->setWriteMCTruthForLowPtHits(m_writeMCTruthForLowPtHits);
      tpcsd->setLowPtCut(m_lowPtCut/dd4hep::MeV*CLHEP::MeV);
      tpcsd->setLowPtMaxHitSeparation(m_lowPtMaxHitSeparation/dd4hep::mm*CLHEP::mm);
      sd = tpcsd;
      info() << "TPC will use TimeProjectionChamberSensitiveDetector" << endmsg;
    }
    else{
      info() << "TPC will use default Geant4TrackerHit SD " << endmsg;
    }
  }

  return sd;
}
