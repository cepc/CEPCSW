#include "DetSimSvc.h"

#include "G4RunManager.hh"

DECLARE_COMPONENT(DetSimSvc)

DetSimSvc::DetSimSvc(const std::string& name, ISvcLocator* svc)
: base_class(name, svc) {
    m_runmgr = nullptr;
}

DetSimSvc::~DetSimSvc() {

}

G4RunManager*
DetSimSvc::getRM() {
    return m_runmgr;
}

StatusCode
DetSimSvc::initializeRM() {
    StatusCode sc;


    G4bool cond = m_runmgr->ConfirmBeamOnCondition();
    if(!cond) {
        return StatusCode::FAILURE;
    }

    m_runmgr->ConstructScoringWorlds();
    m_runmgr->RunInitialization();

    return sc;
}

StatusCode
DetSimSvc::simulateEvent(int i_event) {
    StatusCode sc;

    m_runmgr->ProcessOneEvent(i_event);

    return sc;
}

StatusCode
DetSimSvc::finalizeRM() {
    StatusCode sc;

    m_runmgr->RunTermination();

    return sc;
}

StatusCode
DetSimSvc::initialize() {
    StatusCode sc;

    m_runmgr = new G4RunManager();

    return sc;
}

StatusCode
DetSimSvc::finalize() {
    StatusCode sc;
    delete m_runmgr;

    return sc;
}
