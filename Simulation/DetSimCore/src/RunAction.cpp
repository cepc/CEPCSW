#include "RunAction.h"

#include "G4Run.hh"

RunAction::RunAction(ToolHandleArray<IAnaElemTool>& anatools) 
    : G4UserRunAction(),
      m_anaelemtools(anatools) {

}

RunAction::~RunAction() {

}

void 
RunAction::BeginOfRunAction(const G4Run* aRun)
{

}

void
RunAction::EndOfRunAction(const G4Run* aRun)
{

}


