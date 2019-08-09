#include "DetSimAlg.h"

#include "G4RunManager.hh"

#include "DetectorConstruction.h"

#include "G4PhysListFactory.hh"

#include "PrimaryGeneratorAction.h"

DECLARE_COMPONENT(DetSimAlg)

DetSimAlg::DetSimAlg(const std::string& name, ISvcLocator* pSvcLocator)
: Algorithm(name, pSvcLocator) {
    i_event = -1;
}

StatusCode
DetSimAlg::initialize() {
    StatusCode sc;

    info() << "Initialize DetSimAlg... " << endmsg;

    m_detsimsvc = service("DetSimSvc");
    if (!m_detsimsvc) {
        error() << "Failed to find DetSimSvc. " << endmsg;
        return StatusCode::FAILURE;
    }

    // Initialize the Run Manager
    G4RunManager* runmgr = m_detsimsvc->getRM();
    if (!runmgr) {
        error() << "Failed to get Run Manager. " << endmsg;
        return StatusCode::FAILURE;
    }

    // Detector Construction
    runmgr->SetUserInitialization(new DetectorConstruction());
    // Physics List
    G4VUserPhysicsList *physicsList = nullptr;
    if (m_physics_lists_name.value() == "CEPC") {

    } else {
        G4PhysListFactory *physListFactory = new G4PhysListFactory();
        physicsList = physListFactory->GetReferencePhysList(m_physics_lists_name.value());
    }
    assert(physicsList);
    runmgr->SetUserInitialization(physicsList);
    // Primary Generator Action
    runmgr->SetUserAction(new PrimaryGeneratorAction());

    // after set up the user initialization and user actions, start the initialization.
    m_detsimsvc->initializeRM();

    return sc;
}

StatusCode
DetSimAlg::execute() {
    StatusCode sc;

    m_detsimsvc->simulateEvent(++i_event);

    return sc;
}

StatusCode
DetSimAlg::finalize() {
    StatusCode sc;
    if (!m_detsimsvc) { 
        return StatusCode::FAILURE;
    }
    m_detsimsvc->finalizeRM();


    return sc;
}



