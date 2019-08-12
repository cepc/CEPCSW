#include "AnExampleDetElemTool.h"

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


DECLARE_COMPONENT(AnExampleDetElemTool)

G4LogicalVolume*
AnExampleDetElemTool::getLV() {

    G4Material* Galactic = G4Material::GetMaterial("Galactic");

    G4VSolid* solidAnExample= new G4Box("sAnExample", m_x.value(), m_y.value(), m_z.value());
    G4LogicalVolume* logicAnExample= new G4LogicalVolume( solidAnExample, Galactic, "lAnExample", 0, 0, 0);

    return logicAnExample;
}

StatusCode
AnExampleDetElemTool::initialize() {
    StatusCode sc;
    return sc;
}

StatusCode
AnExampleDetElemTool::finalize() {
    StatusCode sc;
    return sc;
}
