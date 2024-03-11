#include "DetSimSD/DDG4SensitiveDetector.h"

#include "DDG4/Geant4Converter.h"
#include "DD4hep/Segmentations.h"

#include "DD4hep/Printout.h"
#include "DD4hep/Detector.h"

// Geant4 include files
#include "G4Step.hh"
#include "G4PVPlacement.hh"

// ROOT include files
#include "TGeoNode.h"

#include "DD4hep/DD4hepUnits.h"
#include "CLHEP/Units/SystemOfUnits.h"

// convert from CLHEP to DD4hep
static const double MM_2_CM = (dd4hep::millimeter/CLHEP::millimeter);
// convert from DD4hep to CLHEP
static const double CM_2_MM = (CLHEP::millimeter/dd4hep::millimeter);

DDG4SensitiveDetector::DDG4SensitiveDetector(const std::string& name, dd4hep::Detector& description)
    : G4VSensitiveDetector(name), m_detDesc(description),
      m_detector(), m_sensitive(), m_readout() {
    m_detector = description.detector(name);
    m_sensitive = description.sensitiveDetector(name);
    m_readout = m_sensitive.readout();

}

void
DDG4SensitiveDetector::Initialize(G4HCofThisEvent* HCE) {

}

G4bool
DDG4SensitiveDetector::ProcessHits(G4Step* step, G4TouchableHistory*) {
    
    return true;
}

void
DDG4SensitiveDetector::EndOfEvent(G4HCofThisEvent* HCE) {

}

long long
DDG4SensitiveDetector::getVolumeID(const G4Step* aStep) {

    dd4hep::sim::Geant4StepHandler step(aStep);
    dd4hep::sim::Geant4VolumeManager volMgr = dd4hep::sim::Geant4Mapping::instance().volumeManager();
    dd4hep::VolumeID id = volMgr.volumeID(step.preTouchable());

    return id;
}

long long
DDG4SensitiveDetector::getCellID(const G4Step* step) {
    dd4hep::sim::Geant4StepHandler h(step);
    dd4hep::sim::Geant4VolumeManager volMgr = dd4hep::sim::Geant4Mapping::instance().volumeManager();
    dd4hep::VolumeID volID  = volMgr.volumeID(h.preTouchable());

    dd4hep::Segmentation        seg    = m_readout.segmentation();
    if ( seg.isValid() )  {
        G4ThreeVector global = 0.5 * ( h.prePosG4()+h.postPosG4());
        G4ThreeVector local  = h.preTouchable()->GetHistory()->GetTopTransform().TransformPoint(global);
        dd4hep::Position loc(local.x()*MM_2_CM, local.y()*MM_2_CM, local.z()*MM_2_CM);
        dd4hep::Position glob(global.x()*MM_2_CM, global.y()*MM_2_CM, global.z()*MM_2_CM);
        dd4hep::VolumeID cID = seg.cellID(loc,glob,volID);
        return cID;
    }
    return volID;
}

dd4hep::Position
DDG4SensitiveDetector::getNominalPosition(const G4Step* step, long long cellID) {
    if(cellID==0) cellID = getCellID(step);
    dd4hep::sim::Geant4StepHandler h(step);
    dd4hep::Segmentation        seg    = m_readout.segmentation();
    dd4hep::Position loc(0,0,0);
    if ( seg.isValid() ) loc = seg.position(cellID);
    G4ThreeVector local(loc.x()*CM_2_MM, loc.y()*CM_2_MM, loc.z()*CM_2_MM);
    G4ThreeVector global = h.preTouchable()->GetHistory()->GetTopTransform().InverseTransformPoint(local);
    return dd4hep::Position(global.x(), global.y(), global.z());
}
