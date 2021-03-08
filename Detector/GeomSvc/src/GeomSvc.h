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

class dd4hep::DetElement;
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
  
  const gear::ZPlanarParametersImpl*  getVXDParameters() override {return m_vxdParameters;};
  const dd4hep::rec::ZPlanarData* getVXDData() override {return m_vxdData;};
  const dd4hep::rec::ConicalSupportData* getBeamPipeData() override {return m_beamPipeData;};

  const std::map<std::string,double>& getDetParameters(std::string name) override;
  double getDetParameter(std::string set_name, std::string par_name) override;
  TMaterial* getMaterial(std::string name);
  
 private:
  StatusCode convertVXD(dd4hep::DetElement& sub);

    Decoder* getDecoder(const std::string& readout_name) override;

private:
  TGeoNode* FindNode(TGeoNode* mother, char* name);
  // DD4hep XML compact file path
  Gaudi::Property<std::string> m_dd4hep_xmls{this, "compact"};
  
  // 
  dd4hep::Detector* m_dd4hep_geo;


  gear::ZPlanarParametersImpl* m_vxdParameters{nullptr};
  dd4hep::rec::ZPlanarData* m_vxdData{nullptr};
  dd4hep::rec::ConicalSupportData* m_beamPipeData{nullptr};

  //gear::GearParametersImpl* m_vxdInfra;
  std::map<std::string, std::map<std::string,double> > m_detParameters;
  std::map<std::string, TMaterial*> m_materials;
  struct helpLayer {
    double distance =0;
    double offset =0;
    double thickness =0;
    double length =0;
    double width =0;
    double radLength =0;
    double z =0;
    double foam_spacer_radLength =0;
  };
};

#endif // GeomSvc_h
