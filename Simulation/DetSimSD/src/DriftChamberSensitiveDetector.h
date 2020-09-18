#ifndef DriftChamberSensitiveDetector_h
#define DriftChamberSensitiveDetector_h

/*
 * DriftChamberSensitiveDetector is used in Drift Chamber with dE/dx simulator.
 *
 * 19 Sept. 2020, Tao Lin <lintao@ihep.ac.cn>
 */

#include "DetSimSD/DDG4SensitiveDetector.h"
#include "DetSimInterface/IDedxSimTool.h"
#include "GaudiKernel/ToolHandle.h"

class DriftChamberSensitiveDetector: public DDG4SensitiveDetector {
public:
    typedef dd4hep::sim::Geant4TrackerHit TrackerHit;
    typedef G4THitsCollection<TrackerHit> TrackerHitCollection;

public:
    DriftChamberSensitiveDetector(const std::string& name, dd4hep::Detector& description);

    bool setDedxSimTool(ToolHandle<IDedxSimTool>);
    
public:
    // Geant4 interface

    virtual void Initialize(G4HCofThisEvent* HCE);
    virtual G4bool ProcessHits(G4Step* step,G4TouchableHistory* history);
    virtual void EndOfEvent(G4HCofThisEvent* HCE);

protected:

    HitCollection* m_hc;

    // this is passed from SensDetTool
    ToolHandle<IDedxSimTool> m_dedx_simtool;

};

#endif
