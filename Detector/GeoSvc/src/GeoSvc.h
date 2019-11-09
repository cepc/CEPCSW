#ifndef GeoSvc_h
#define GeoSvc_h

// Interface
#include "DetInterface/IGeoSvc.h"

// Gaudi
#include "GaudiKernel/IIncidentListener.h"
#include "GaudiKernel/IIncidentSvc.h"
#include "GaudiKernel/Incident.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/Service.h"
#include "GaudiKernel/ServiceHandle.h"

// DD4Hep
#include "DD4hep/Detector.h"

class GeoSvc: public extends<Service, IGeoSvc> {
public:
    GeoSvc(const std::string& name, ISvcLocator* svc);
    ~GeoSvc();

    // Service
    StatusCode initialize() override;
    StatusCode finalize() override;

    // IGeoSvc
    dd4hep::DetElement getDD4HepGeo() override;
    dd4hep::Detector* lcdd() override;

private:

    // DD4hep XML compact file path
    Gaudi::Property<std::string> m_dd4hep_xmls{this, "compact"};

    // 
    dd4hep::Detector* m_dd4hep_geo;
};


#endif GeoSvc_h
