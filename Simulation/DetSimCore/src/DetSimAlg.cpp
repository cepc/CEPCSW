#include "DetSimAlg.h"
#include "GaudiKernel/IEventProcessor.h"
#include "GaudiKernel/IAppMgrUI.h"
#include "GaudiKernel/GaudiException.h"

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

#include "DetectorConstruction.h"
#include "G4PhysListFactory.hh"
#include "PrimaryGeneratorAction.h"

#include "ActionInitialization.h"

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

    // User Actions
    runmgr->SetUserInitialization(new ActionInitialization());

    // Vis Mac
    bool hasVis = false;
    G4VisManager* visManager = nullptr;
    G4UIExecutive* ui = nullptr;
    if (not m_vis_macs.value().empty()) {
        hasVis = true;
        info() << "Start Geant4 Visualization." << endmsg;
        // initialize Vis
        visManager = new G4VisExecutive;
        char* argv[1] = {strdup("geant4Vis")};
        ui = new G4UIExecutive(1,argv);
    }
    for (auto vismac: m_vis_macs.value()) {
        G4UImanager *UImanager = G4UImanager::GetUIpointer();
        std::string command = "/control/execute ";
        UImanager->ApplyCommand( command + vismac);
    }

    // Run Mac & Run Cmds
    for (auto runmac: m_run_macs.value()) {
        G4UImanager *UImanager = G4UImanager::GetUIpointer();
        std::string command = "/control/execute ";
        UImanager->ApplyCommand( command + runmac);
    }
    for (auto runcmd: m_run_cmds.value()) {
        G4UImanager *UImanager = G4UImanager::GetUIpointer();
        UImanager->ApplyCommand(runcmd);
    }

    // if has vis, we stop here
    if (hasVis) {

        ui->SessionStart();

        delete ui;
        delete visManager;

        // fixme: how to stop the run?
        // auto ep = serviceLocator()->as<IEventProcessor>();
        // if ( !ep ) {
        //     error() << "Cannot get IEventProcessor" << endmsg;
        //     return StatusCode::FAILURE;
        // }
        // ep->stopRun();

        // auto ep = serviceLocator()->as<IAppMgrUI>();
        // if ( !ep ) {
        //     error() << "Cannot get IAppMgrUI" << endmsg;
        //     return StatusCode::FAILURE;
        // }
        // ep->terminate();
        // throw GaudiException( "Stopping loop", "Stop", StatusCode::SUCCESS);
        return StatusCode::FAILURE;
    }

    // Initialize G4 Kernel
    runmgr->Initialize();

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



