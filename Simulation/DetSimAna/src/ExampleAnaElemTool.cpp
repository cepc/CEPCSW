#include "ExampleAnaElemTool.h"

#include "G4Event.hh"

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
ExampleAnaElemTool::EndOfEventAction(const G4Event*) {

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


