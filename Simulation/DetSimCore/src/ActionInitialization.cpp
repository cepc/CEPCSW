#include "ActionInitialization.h"

#include "RunAction.h"

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

}
