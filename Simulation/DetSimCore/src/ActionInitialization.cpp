#include "ActionInitialization.h"

#include "RunAction.h"
#include "EventAction.h"
#include "TrackingAction.h"
#include "SteppingAction.h"

ActionInitialization::ActionInitialization(ToolHandleArray<IAnaElemTool>& anatools)
    : G4VUserActionInitialization(),
      m_anaelemtools(anatools) {

}

ActionInitialization::~ActionInitialization() {

}

void
ActionInitialization::BuildForMaster() const {

}

void
ActionInitialization::Build() const {


    RunAction* runAction = new RunAction(m_anaelemtools);
    SetUserAction(runAction);

    EventAction* eventAction = new EventAction(m_anaelemtools);
    SetUserAction(eventAction);

    TrackingAction* trackingAction = new TrackingAction(m_anaelemtools);
    SetUserAction(trackingAction);

    SteppingAction* steppingAction = new SteppingAction(m_anaelemtools);
    SetUserAction(steppingAction);

}
