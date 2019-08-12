#ifndef DetectorConstruction_h
#define DetectorConstruction_h

#include "GaudiKernel/ToolHandle.h"
#include "DetSimInterface/IDetElemTool.h"

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4OpticalSurface.hh"
#include "G4Material.hh"

// A concrete detector construction class.
// The base class is Geant4's G4VUserDetectorConstruction only.
// Another Gaudi tool is used to configure & create this object.

class DetectorConstruction: public G4VUserDetectorConstruction {

public:
    DetectorConstruction(ToolHandle<IDetElemTool>& root_elem);
    ~DetectorConstruction();
public:
    G4VPhysicalVolume* Construct();

private:
    ToolHandle<IDetElemTool>& m_root_detelem;
};

#endif
