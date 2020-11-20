#ifndef DDG4SensitiveDetector_h
#define DDG4SensitiveDetector_h

/*
 * In order to access ID from DDG4, some utilities are necessary to retrieve information
 * from DDG4. This base class defines such interfaces and utilities.
 * 
 * Refer to the class DDG4/include/DDG4/Geant4SensitiveDetector.h for some APIs usage.
 *
 * We keep to reuse some types already defined in DDG4:
 * - Geant4Hits
 *
 * -- 12 June 2020, Tao Lin <lintao@ihep.ac.cn>
 */


#include "DD4hep/Detector.h"
#include "DDG4/Geant4Hits.h"

#include "G4Step.hh"
#include "G4HCofThisEvent.hh"
#include "G4TouchableHistory.hh"
#include "G4VSensitiveDetector.hh"
#include "G4THitsCollection.hh"


class DDG4SensitiveDetector: public G4VSensitiveDetector {
public:
    typedef dd4hep::sim::Geant4Hit                    Geant4Hit;
    typedef G4THitsCollection<dd4hep::sim::Geant4Hit> HitCollection;
    typedef dd4hep::sim::Geant4Hit::Contribution      HitContribution;

public:
    DDG4SensitiveDetector(const std::string& name, dd4hep::Detector& description);

public:
    // Geant4 interface

    virtual void Initialize(G4HCofThisEvent* HCE);
    virtual G4bool ProcessHits(G4Step* step,G4TouchableHistory* history);
    virtual void EndOfEvent(G4HCofThisEvent* HCE);

public:
    // DDG4 utilities
    /// Returns the volumeID of the sensitive volume corresponding to the step -
    /// combining the VolIDS of the complete geometry path (Geant4TouchableHistory)
    //  from the current sensitive volume to the world volume
    virtual long long getVolumeID(const G4Step* step);

    /// Returns the volumeID of the sensitive volume corresponding to the step -
    /// combining the VolIDS of the complete geometry path (Geant4TouchableHistory)
    //  from the current sensitive volume to the world volume
    virtual long long getCellID(const G4Step* step);


protected:
    /// Reference to the detector description object
    dd4hep::Detector& m_detDesc;

    /// Reference to the detector element describing this sensitive element
    dd4hep::DetElement m_detector;

    /// Reference to the sensitive detector element
    dd4hep::SensitiveDetector m_sensitive;

    /// Reference to the readout structure
    dd4hep::Readout m_readout;

};

#endif
