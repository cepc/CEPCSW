#include "DetectorConstruction.h"

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

#include "G4Isotope.hh"
#include "G4Element.hh"
#include "G4MaterialPropertiesTable.hh"

#include "G4Box.hh"
#include "G4Sphere.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SDManager.hh"
#include "G4Region.hh"
#include "G4RegionStore.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4OpticalSurface.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4ios.hh"


DetectorConstruction::DetectorConstruction(ToolHandle<IDetElemTool>& root_elem) 
    : m_root_detelem(root_elem) {

}

DetectorConstruction::~DetectorConstruction() {

}

G4VPhysicalVolume* 
DetectorConstruction::Construct() {
    // =======================================================================
    // Materials
    // =======================================================================
    bool any_warnings = false;

    G4Material* Galactic = G4Material::GetMaterial("Galactic", any_warnings);
    if (not Galactic) { 
        Galactic = new G4Material("Galactic", 1., 1.01*g/mole, universe_mean_density, kStateGas, 2.73*kelvin, 3.e-18*pascal);
    }


    // =======================================================================
    // World
    // =======================================================================
    G4LogicalVolume* logicWorld = m_root_detelem->getLV();

    G4VPhysicalVolume* physiWorld = new G4PVPlacement(0,               // no rotation
                                                      G4ThreeVector(), // at (0,0,0)
                                                      logicWorld,      // its logical volume
                                                      "pWorld",         // its name
                                                      0,               // its mother  volume
                                                      false,           // no boolean operations
                                                      0);              // no field specific to volume

    return physiWorld;

}

void
DetectorConstruction::ConstructSDandField() {

}
