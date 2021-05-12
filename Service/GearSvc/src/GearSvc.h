#ifndef GEAR_SVC_H
#define GEAR_SVC_H

#include "GearSvc/IGearSvc.h"
#include <GaudiKernel/Service.h>
#include "DD4hep/Detector.h"
class dd4hep::DetElement;
class TGeoNode;

class GearSvc : public extends<Service, IGearSvc>
{
    public:
        GearSvc(const std::string& name, ISvcLocator* svc);
        ~GearSvc();

        gear::GearMgr* getGearMgr() override;

        StatusCode initialize() override;
        StatusCode finalize() override;

    private:
	StatusCode convertBeamPipe(dd4hep::DetElement& pipe);
	StatusCode convertVXD(dd4hep::DetElement& vxd);
	StatusCode convertSIT(dd4hep::DetElement& sit);
	StatusCode convertTPC(dd4hep::DetElement& tpc);
	StatusCode convertDC (dd4hep::DetElement& dc);
	StatusCode convertSET(dd4hep::DetElement& set);
	StatusCode convertFTD(dd4hep::DetElement& ftd);
	TGeoNode* FindNode(TGeoNode* mother, char* name);

        Gaudi::Property<std::string> m_gearFile{this, "GearXMLFile", ""};

        gear::GearMgr* m_gearMgr;
};

#endif
