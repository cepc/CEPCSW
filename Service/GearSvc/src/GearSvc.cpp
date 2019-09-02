#include "GearSvc.h"
#include "gearxml/GearXML.h"
#include "gearimpl/GearMgrImpl.h"

DECLARE_COMPONENT(GearSvc)

GearSvc::GearSvc(const std::string& name, ISvcLocator* svc)
    : base_class(name, svc),
      m_gearMgr(nullptr)
{
}

GearSvc::~GearSvc()
{
}

gear::GearMgr* GearSvc::getGearMgr()
{
    return m_gearMgr;
}

StatusCode GearSvc::initialize()
{
    if ( m_gearFile.size() > 0 ) {
        info() << "instantiated GEAR from file " << m_gearFile << endmsg;
        m_gearMgr = gear::GearXML(m_gearFile).createGearMgr();
    }
    else {
        warning() << "no GEAR XML file given ..." << endmsg;
        m_gearMgr = new gear::GearMgrImpl;
    }

    return StatusCode::SUCCESS;
}

StatusCode GearSvc::finalize()
{
    if ( m_gearMgr ) {
        delete m_gearMgr;
        m_gearMgr = nullptr;
    }

    return StatusCode::SUCCESS;
}
