#include "ExampleAnaElemTool.h"

#include "G4Event.hh"

#include "DD4hep/Detector.h"
#include "DD4hep/Plugins.h"
#include "DDG4/Geant4Converter.h"
#include "DDG4/Geant4Mapping.h"


DECLARE_COMPONENT(ExampleAnaElemTool)

void
ExampleAnaElemTool::BeginOfRunAction(const G4Run*) {
    G4cout << "Begin Run of detector simultion..." << G4endl;
}

void
ExampleAnaElemTool::EndOfRunAction(const G4Run*) {
    G4cout << "End Run of detector simultion..." << G4endl;
}

void
ExampleAnaElemTool::BeginOfEventAction(const G4Event* anEvent) {
    msg() << "Event " << anEvent->GetEventID() << endmsg;
}

void
ExampleAnaElemTool::EndOfEventAction(const G4Event* anEvent) {

    // save all data

    // readout defined in DD4hep
    auto lcdd = &(dd4hep::Detector::getInstance());
    auto allReadouts = lcdd->readouts();

    for (auto& readout : allReadouts) {
        info() << "Readout " << readout.first << endmsg;
    }

    // retrieve the hit collections
    G4HCofThisEvent* collections = anEvent->GetHCofThisEvent();
    if (!collections) {
        warning() << "No collections found. " << endmsg;
        return;
    }
    int Ncol = collections->GetNumberOfCollections();
    for (int icol = 0; icol < Ncol; ++icol) {
        G4VHitsCollection* collect = collections->GetHC(icol);
        if (!collect) {
            warning() << "Collection iCol " << icol << " is missing" << endmsg;
            continue;
        }
        info() << "Collection " << collect->GetName() << endmsg;
    }
}

void
ExampleAnaElemTool::PreUserTrackingAction(const G4Track*) {

}

void
ExampleAnaElemTool::PostUserTrackingAction(const G4Track*) {

}

void
ExampleAnaElemTool::UserSteppingAction(const G4Step*) {

}

StatusCode
ExampleAnaElemTool::initialize() {
    StatusCode sc;

    return sc;
}

StatusCode
ExampleAnaElemTool::finalize() {
    StatusCode sc;

    return sc;
}


