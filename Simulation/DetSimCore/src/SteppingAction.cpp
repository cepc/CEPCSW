#include "SteppingAction.h"

SteppingAction::SteppingAction(ToolHandleArray<IAnaElemTool>& anatools)
    : m_anaelemtools(anatools) {

}

SteppingAction::~SteppingAction() {

}

void
SteppingAction::UserSteppingAction(const G4Step* aStep) {

}

