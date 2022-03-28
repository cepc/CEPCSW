// *********************************************************
//
// $Id: GenericTrackerSensitiveDetector.hh,v 1.0 2022/03/27

#ifndef GenericTrackerSensitiveDetector_h
#define GenericTrackerSensitiveDetector_h

#include "DetSimSD/DDG4SensitiveDetector.h"
#include "DDG4/Defs.h"

class GenericTrackerSensitiveDetector: public DDG4SensitiveDetector {
 public:
  GenericTrackerSensitiveDetector(const std::string& name, dd4hep::Detector& description);
  
  void Initialize(G4HCofThisEvent* HCE);
  G4bool ProcessHits(G4Step* step, G4TouchableHistory* history);
  void EndOfEvent(G4HCofThisEvent* HCE);
  
 protected:

  HitCollection* m_hc = nullptr;
  
};
#endif
