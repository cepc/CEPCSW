#include "GeomSvc.h"
#include "gearimpl/GearParametersImpl.h"
#include "TMath.h"
#include "TMaterial.h"
#include "CLHEP/Units/SystemOfUnits.h"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4NistManager.hh"

#include "DD4hep/Detector.h"
#include "DD4hep/Plugins.h"
#include "DDG4/Geant4Converter.h"
#include "DDG4/Geant4Mapping.h"
#include "DDRec/DetectorData.h"

#include <iomanip>
#include <iostream>

DECLARE_COMPONENT(GeomSvc)

GeomSvc::GeomSvc(const std::string& name, ISvcLocator* svc)
: base_class(name, svc), m_dd4hep_geo(nullptr){

}

GeomSvc::~GeomSvc() {

}

StatusCode
GeomSvc::initialize() {
  StatusCode sc = Service::initialize();

  m_dd4hep_geo = &(dd4hep::Detector::getInstance());
  // if failed to load the compact, a runtime error will be thrown.
  m_dd4hep_geo->fromCompact(m_dd4hep_xmls.value());
  
  return sc;
}

StatusCode
GeomSvc::finalize() {
  StatusCode sc;

  dd4hep::Detector::destroyInstance();

  return sc;
}

dd4hep::DetElement
GeomSvc::getDD4HepGeo() {
    if (lcdd()) {
        return lcdd()->world();
    }
    return dd4hep::DetElement();
}

dd4hep::Detector*
GeomSvc::lcdd() {
    return m_dd4hep_geo;
}


IGeomSvc::Decoder*
GeomSvc::getDecoder(const std::string& readout_name) {

    IGeomSvc::Decoder* decoder = nullptr;

    if (!lcdd()) {
        error() << "Failed to get lcdd()" << endmsg;
        return decoder;
    }

    auto readouts = m_dd4hep_geo->readouts();
    if (readouts.find(readout_name) == readouts.end()) {
        error() << "Failed to find readout name '" << readout_name << "'"
                << " in DD4hep::readouts. "
                << endmsg;
        return decoder;
    }
    
    dd4hep::Readout readout = lcdd()->readout(readout_name);
    auto m_idspec = readout.idSpec(); 

    decoder = m_idspec.decoder();

    if (!decoder) {
        error() << "Failed to get the decoder with readout '"
                << readout_name << "'" << endmsg;
    }

    return decoder;

}
