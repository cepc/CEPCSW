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

#include "DD4hep/Detector.h"
#include "DDG4/Geant4Converter.h"
#include "DDG4/Geant4Mapping.h"

DECLARE_COMPONENT(AnExampleDetElemTool)

G4LogicalVolume*
AnExampleDetElemTool::getLV() {

    G4Material* Galactic = G4Material::GetMaterial("Galactic");

    G4VSolid* solidAnExample= new G4Box("sAnExample", m_x.value(), m_y.value(), m_z.value());
    G4LogicalVolume* logicAnExample= new G4LogicalVolume( solidAnExample, Galactic, "lAnExample", 0, 0, 0);

    // Following is an example to get the DD4hep volume
    dd4hep::Detector* dd4hep_geo = &(dd4hep::Detector::getInstance());
    dd4hep_geo->fromCompact(m_dd4hep_xmls.value());
    dd4hep::DetElement world = dd4hep_geo->world();
    dd4hep::sim::Geant4Converter conv((*dd4hep_geo), dd4hep::DEBUG);

    dd4hep::sim::Geant4GeometryInfo* geo_info = conv.create(world).detach();
    dd4hep::sim::Geant4Mapping&  g4map = dd4hep::sim::Geant4Mapping::instance();
    g4map.attach(geo_info);
    // All volumes are deleted in ~G4PhysicalVolumeStore()
    G4VPhysicalVolume* m_world = geo_info->world();
    G4LogicalVolume* logicDD4hepExample = m_world->GetLogicalVolume();

    if (logicDD4hepExample) {
        new G4PVPlacement(0,                   // no rotation
                          G4ThreeVector(),     // at (0,0,0)
                          logicDD4hepExample,  // logical volume
                          "lDD4hepExampleDetElem", // name
                          logicAnExample,      // mother volume
                          false,               // no boolean operations
                          0);                  // no field
    } else {
        warning() << "Can't Find the logical volume lDD4hepExampleDetElem " << std::endl;
    }


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
