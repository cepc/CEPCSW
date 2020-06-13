#ifndef CaloSensitiveDetector_h
#define CaloSensitiveDetector_h

/*
 * This is an implementation of Calo SD.
 *
 * -- 13 June 2020, Tao Lin <lintao@ihep.ac.cn>
 */

#include "DetSimSD/DDG4SensitiveDetector.h"

class CaloSensitiveDetector: public DDG4SensitiveDetector {
public:
    typedef dd4hep::sim::Geant4CalorimeterHit CalorimeterHit;
    typedef G4THitsCollection<CalorimeterHit> CaloHitCollection;

public:
    CaloSensitiveDetector(const std::string& name, dd4hep::Detector& description);

public:
    // Geant4 interface

    virtual void Initialize(G4HCofThisEvent* HCE);
    virtual G4bool ProcessHits(G4Step* step,G4TouchableHistory* history);
    virtual void EndOfEvent(G4HCofThisEvent* HCE);

protected:
    CalorimeterHit* find(const HitCollection*, const dd4hep::sim::HitCompare<CalorimeterHit>&);
    
protected:

    HitCollection* m_hc;
};


#endif
