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
