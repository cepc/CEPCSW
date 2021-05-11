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

// Field
#include "G4UniformMagField.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"

#include "DD4hep/Detector.h"
#include "DD4hep/Plugins.h"
#include "DDG4/Geant4Converter.h"
#include "DDG4/Geant4Mapping.h"
#include "DDG4/Geant4Field.h"

DECLARE_COMPONENT(AnExampleDetElemTool)

G4LogicalVolume*
AnExampleDetElemTool::getLV() {

    // G4Material* Galactic = G4Material::GetMaterial("Galactic");

    // G4VSolid* solidAnExample= new G4Box("sAnExample", m_x.value(), m_y.value(), m_z.value());
    // G4LogicalVolume* logicAnExample= new G4LogicalVolume( solidAnExample, Galactic, "lAnExample", 0, 0, 0);

    // Following is an example to get the DD4hep volume
    // dd4hep::Detector* dd4hep_geo = &(dd4hep::Detector::getInstance());
    // dd4hep_geo->fromCompact(m_dd4hep_xmls.value());
    // dd4hep::DetElement world = dd4hep_geo->world();

    dd4hep::Detector* dd4hep_geo = m_geosvc->lcdd();
    dd4hep::DetElement world = m_geosvc->getDD4HepGeo();

    dd4hep::sim::Geant4Converter conv((*dd4hep_geo), dd4hep::DEBUG);

    dd4hep::sim::Geant4GeometryInfo* geo_info = conv.create(world).detach();
    dd4hep::sim::Geant4Mapping&  g4map = dd4hep::sim::Geant4Mapping::instance();
    g4map.attach(geo_info);
    // All volumes are deleted in ~G4PhysicalVolumeStore()
    G4VPhysicalVolume* m_world = geo_info->world();
    G4LogicalVolume* logicDD4hepExample = m_world->GetLogicalVolume();

    if (logicDD4hepExample) {
        // new G4PVPlacement(0,                   // no rotation
        //                   G4ThreeVector(),     // at (0,0,0)
        //                   logicDD4hepExample,  // logical volume
        //                   "lDD4hepExampleDetElem", // name
        //                   logicAnExample,      // mother volume
        //                   false,               // no boolean operations
        //                   0);                  // no field
    } else {
        warning() << "Can't Find the logical volume lDD4hepExampleDetElem " << std::endl;
    }


    // return logicAnExample;
    return logicDD4hepExample;
}

void
AnExampleDetElemTool::ConstructSDandField() {
    //
    // Construct SD using DD4hep.
    // Refer to FCCSW/Detector/DetComponents/src/
    // 

    typedef std::set<const TGeoVolume*> VolSet;
    typedef std::map<dd4hep::SensitiveDetector, VolSet> _SV;
    dd4hep::sim::Geant4GeometryInfo* p = dd4hep::sim::Geant4Mapping::instance().ptr();
    _SV& vols = p->sensitives;

    auto lcdd = m_geosvc->lcdd();

    for (_SV::const_iterator iv = vols.begin(); iv != vols.end(); ++iv) {
        dd4hep::SensitiveDetector sd = (*iv).first;
        std::string typ = sd.type(), nam = sd.name();

        info() << "Type/Name: "
               << typ << "/" << nam
               << endmsg;
        // continue;
        // Sensitive detectors are deleted in ~G4SDManager
        G4VSensitiveDetector* g4sd = nullptr;

        // try to use SD tool to find the SD
        if (!g4sd) {
            if (typ=="calorimeter") {
                m_calo_sdtool = ToolHandle<ISensDetTool>("CalorimeterSensDetTool");
                if (m_calo_sdtool) {
                    info() << "Find the CalorimeterSensDetTool." << endmsg;
                    g4sd = m_calo_sdtool->createSD(nam);
                    info() << "create g4SD: " << g4sd << endmsg;
                }
            } else if (typ=="tracker") {

                // if drift chamber
                if (nam == "DriftChamber") {
                    m_driftchamber_sdtool = ToolHandle<ISensDetTool>("DriftChamberSensDetTool");
                    if (m_driftchamber_sdtool) {
                        info() << "Find the DriftChamberSensDetTool" << endmsg;
                        g4sd = m_driftchamber_sdtool->createSD(nam);
                    } else {
                        warning() << "DriftChamberSensDetTool is not found. " << endmsg;
                    }
                }
		else if (nam == "TPC") {
		  m_tpc_sdtool = ToolHandle<ISensDetTool>("TimeProjectionChamberSensDetTool");
		  if (m_tpc_sdtool) {
		    info() << "Find the TimeProjectionChamberSensDetTool" << endmsg;
		    g4sd = m_tpc_sdtool->createSD(nam);
		  }
		  else {
		    warning() << "TimeProjectionChamberSensDetTool is not found, and default tracker SD will be used" << endmsg;
		  }
		}

            }
        }
        
        if (!g4sd) {
            g4sd = dd4hep::PluginService::Create<G4VSensitiveDetector*>(typ, nam, lcdd);
        }

        if (g4sd == nullptr) {
            std::string tmp = typ;
            tmp[0] = ::toupper(tmp[0]);
            typ = "Geant4" + tmp;
            g4sd = dd4hep::PluginService::Create<G4VSensitiveDetector*>(typ, nam, lcdd);
            if (g4sd == nullptr) {
                dd4hep::PluginDebug dbg;
                g4sd = dd4hep::PluginService::Create<G4VSensitiveDetector*>(typ, nam, lcdd);
                if (g4sd == nullptr) {
                    throw std::runtime_error("ConstructSDandField: FATAL Failed to "
                                             "create Geant4 sensitive detector " +
                                             nam + " of type " + typ + ".");
                }
            }
        }
        g4sd->Activate(true);
        G4SDManager::GetSDMpointer()->AddNewDetector(g4sd);
        const VolSet& sens_vols = (*iv).second;
        for (VolSet::const_iterator i = sens_vols.begin(); i != sens_vols.end(); ++i) {
            const TGeoVolume* vol = *i;
            G4LogicalVolume* g4v = p->g4Volumes[vol];
            if (g4v == nullptr) {
                throw std::runtime_error("ConstructSDandField: Failed to access G4LogicalVolume for SD " + nam + " of type " +
                                         typ + ".");
            }
            info() << " -> Adding " << g4v->GetName() << endmsg;
            G4SDManager::GetSDMpointer()->AddNewDetector(g4sd);
            g4v->SetSensitiveDetector(g4sd);
        }
    }

    // =======================================================================
    // Construct Field
    // =======================================================================
    // TODO: integrate the field between DD4hep and Geant4
    // Note:
    //   DD4hep provides the parameters of fields
    //   Geant4 will setup the field based on the DD4hep fields
    
    // Related Examples:
    // - G4: G4GlobalMagFieldMessenger.cc

    G4FieldManager* fieldManager
        = G4TransportationManager::GetTransportationManager()->GetFieldManager();

    // // Below is a uniform B-field
    // G4ThreeVector value(0,0,3.*tesla);
    // G4UniformMagField* mag_field = new G4UniformMagField(value);

    // DDG4 based B-field
    dd4hep::OverlayedField fld  = lcdd->field();
    G4MagneticField* mag_field  = new dd4hep::sim::Geant4Field(fld);

    fieldManager->SetDetectorField(mag_field);
    fieldManager->CreateChordFinder(mag_field);


}

StatusCode
AnExampleDetElemTool::initialize() {
    StatusCode sc;

    m_geosvc = service<IGeomSvc>("GeomSvc");
    if (!m_geosvc) {
        error() << "Failed to find GeomSvc." << endmsg;
        return StatusCode::FAILURE;
    }

    m_calo_sdtool = ToolHandle<ISensDetTool>("CalorimeterSensDetTool");
    m_driftchamber_sdtool = ToolHandle<ISensDetTool>("DriftChamberSensDetTool");
    m_tpc_sdtool = ToolHandle<ISensDetTool>("TimeProjectionChamberSensDetTool");

    return sc;
}

StatusCode
AnExampleDetElemTool::finalize() {
    StatusCode sc;
    return sc;
}
