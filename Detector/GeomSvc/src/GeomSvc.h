#ifndef GeomSvc_h
#define GeomSvc_h

// Interface
#include "DetInterface/IGeomSvc.h"

// Gaudi
#include "GaudiKernel/IIncidentListener.h"
#include "GaudiKernel/IIncidentSvc.h"
#include "GaudiKernel/Incident.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/Service.h"
#include "GaudiKernel/ServiceHandle.h"

// DD4Hep
#include "DD4hep/Detector.h"

#include <gear/GEAR.h>
#include <gearimpl/ZPlanarParametersImpl.h>
#include <gearimpl/GearParametersImpl.h>

class TGeoNode;

class GeomSvc: public extends<Service, IGeomSvc> {
 public:
  GeomSvc(const std::string& name, ISvcLocator* svc);
  ~GeomSvc();
  
  // Service
  StatusCode initialize() override;
  StatusCode finalize() override;
  
  // IGeomSvc
  dd4hep::DetElement getDD4HepGeo() override;
  dd4hep::Detector* lcdd() override;
  
 private:
  Decoder* getDecoder(const std::string& readout_name) override;
    
private:
  // DD4hep XML compact file path
  Gaudi::Property<std::string> m_dd4hep_xmls{this, "compact"};
  
  // 
  dd4hep::Detector* m_dd4hep_geo;
};

#endif // GeomSvc_h
