#include "ActionInitialization.h"

#include "RunAction.h"
#include "EventAction.h"
#include "TrackingAction.h"
#include "SteppingAction.h"

ActionInitialization::ActionInitialization()
    : G4VUserActionInitialization() {

}

ActionInitialization::~ActionInitialization() {

}

void
ActionInitialization::BuildForMaster() const {

}

void
ActionInitialization::Build() const {


    RunAction* runAction = new RunAction();
    SetUserAction(runAction);

    EventAction* eventAction = new EventAction();
    SetUserAction(eventAction);

    TrackingAction* trackingAction = new TrackingAction();
    SetUserAction(trackingAction);

    SteppingAction* steppingAction = new SteppingAction();
    SetUserAction(steppingAction);

}
