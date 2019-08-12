#include "WorldDetElemTool.h"

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


DECLARE_COMPONENT(WorldDetElemTool)

G4LogicalVolume*
WorldDetElemTool::getLV() {

    G4Material* Galactic = G4Material::GetMaterial("Galactic");

    G4VSolid* solidWorld= new G4Box("sWorld", 60*m, 60*m, 60*m);
    G4LogicalVolume* logicWorld= new G4LogicalVolume( solidWorld, Galactic, "lWorld", 0, 0, 0);

    return logicWorld;
}

StatusCode
WorldDetElemTool::initialize() {
    StatusCode sc;
    return sc;
}

StatusCode
WorldDetElemTool::finalize() {
    StatusCode sc;
    return sc;
}
