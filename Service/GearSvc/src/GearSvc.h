#ifndef GEAR_SVC_H
#define GEAR_SVC_H

#include "GearSvc/IGearSvc.h"
#include <GaudiKernel/Service.h>

class GearSvc : public extends<Service, IGearSvc>
{
    public:
        GearSvc(const std::string& name, ISvcLocator* svc);
        ~GearSvc();

        gear::GearMgr* getGearMgr() override;

        StatusCode initialize() override;
        StatusCode finalize() override;

    private:

        Gaudi::Property<std::string> m_gearFile{this, "GearXMLFile", ""};

        gear::GearMgr* m_gearMgr;
};

#endif
