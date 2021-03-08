#include "DetSimAlg.h"
#include "GaudiKernel/IEventProcessor.h"
#include "GaudiKernel/IAppMgrUI.h"
#include "GaudiKernel/GaudiException.h"
#include "GaudiKernel/IRndmEngine.h"

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

#include "Randomize.hh"
#include "DetectorConstruction.h"
#include "G4PhysListFactory.hh"
#include "G4EmParameters.hh"
#include "G4StepLimiterPhysics.hh"
#include "G4FastSimulationPhysics.hh"
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

    // Initialize random seed
    if (not m_randomSeeds.empty()) {
        randSvc()->engine()->setSeeds( m_randomSeeds );
    }

    info() << "Random Seed is initialized to "
           << G4Random::getTheSeed()
           << " in Geant4" << endmsg;

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
    m_root_detelem = ToolHandle<IDetElemTool>(m_root_det_elem.value());

    for (auto fastsimname: m_fast_simnames) {
        info() << "Fast Sim Tool: " << fastsimname << endmsg;
        m_fast_simtools.push_back(fastsimname);
    }

    runmgr->SetUserInitialization(new DetectorConstruction(m_root_detelem, m_fast_simtools));
    // Physics List
    G4VUserPhysicsList *physicsList = nullptr;
    if (m_physics_lists_name.value() == "CEPC") {

    } else {
        G4PhysListFactory *physListFactory = new G4PhysListFactory();
        G4VModularPhysicsList* modularPhysicsList = physListFactory->GetReferencePhysList(m_physics_lists_name.value());

        // PAI model
        G4EmParameters::Instance()->AddPAIModel("all","DriftChamberRegion","pai");
        // G4EmParameters::Instance()->AddPAIModel("all","DriftChamberRegion","pai_photon");

        // register addition physics list
        modularPhysicsList->RegisterPhysics(new G4StepLimiterPhysics());

        // register fastsim physics
        G4FastSimulationPhysics* fastsim_physics = new G4FastSimulationPhysics();
        fastsim_physics->BeVerbose();
        fastsim_physics->ActivateFastSimulation("e-");
        fastsim_physics->ActivateFastSimulation("e+");
        fastsim_physics->ActivateFastSimulation("gamma");
        fastsim_physics->ActivateFastSimulation("mu-");
        fastsim_physics->ActivateFastSimulation("pi0");
        fastsim_physics->ActivateFastSimulation("pi+");
        fastsim_physics->ActivateFastSimulation("pi-");
        fastsim_physics->ActivateFastSimulation("tau+");
        fastsim_physics->ActivateFastSimulation("tau-");
        fastsim_physics->ActivateFastSimulation("K0");
        fastsim_physics->ActivateFastSimulation("K-");
        fastsim_physics->ActivateFastSimulation("K+");
        fastsim_physics->ActivateFastSimulation("Z0");
        fastsim_physics->ActivateFastSimulation("W-");
        fastsim_physics->ActivateFastSimulation("h0");
        fastsim_physics->ActivateFastSimulation("nu_mu");
        fastsim_physics->ActivateFastSimulation("nu_ebar");
        modularPhysicsList->RegisterPhysics(fastsim_physics);

        physicsList = modularPhysicsList;
    }
    assert(physicsList);
    runmgr->SetUserInitialization(physicsList);
    // Primary Generator Action
    if (!m_prim_cnvtool) {
        error() << "Failed to get the primary cnvtool." << endmsg;
        return StatusCode::FAILURE;
    }
    runmgr->SetUserAction(new PrimaryGeneratorAction(m_prim_cnvtool));

    // User Actions
    for (auto anaelem: m_ana_elems.value()) {
        m_anaelemtools.push_back(anaelem);
    }
    runmgr->SetUserInitialization(new ActionInitialization(m_anaelemtools));

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



