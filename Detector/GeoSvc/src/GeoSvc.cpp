#include "GeoSvc.h"

#include "DD4hep/Detector.h"
#include "DD4hep/Plugins.h"
#include "DDG4/Geant4Converter.h"
#include "DDG4/Geant4Mapping.h"

DECLARE_COMPONENT(GeoSvc)

GeoSvc::GeoSvc(const std::string& name, ISvcLocator* svc)
: base_class(name, svc), m_dd4hep_geo(nullptr) {

}

GeoSvc::~GeoSvc() {

}

StatusCode
GeoSvc::initialize() {
    StatusCode sc = Service::initialize();

    m_dd4hep_geo = &(dd4hep::Detector::getInstance());
    // if failed to load the compact, a runtime error will be thrown.
    m_dd4hep_geo->fromCompact(m_dd4hep_xmls.value());

    return sc;
}

StatusCode
GeoSvc::finalize() {
    StatusCode sc;

    return sc;
}

dd4hep::DetElement
GeoSvc::getDD4HepGeo() {
    if (lcdd()) {
        return lcdd()->world();
    }
    return dd4hep::DetElement();
}

dd4hep::Detector*
GeoSvc::lcdd() {
    return m_dd4hep_geo;
}

IGeoSvc::Decoder*
GeoSvc::getDecoder(const std::string& readout_name) {

    IGeoSvc::Decoder* decoder = nullptr;

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
