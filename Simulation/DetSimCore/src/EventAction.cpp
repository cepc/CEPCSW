#include "EventAction.h"

EventAction::EventAction(ToolHandleArray<IAnaElemTool>& anatools) 
    : G4UserEventAction(),
      m_anaelemtools(anatools) {

}

EventAction::~EventAction() {

}

void
EventAction::BeginOfEventAction(const G4Event*) {

}

void
EventAction::EndOfEventAction(const G4Event*) {

}

