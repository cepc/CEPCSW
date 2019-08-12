#include "TrackingAction.h"

TrackingAction::TrackingAction(ToolHandleArray<IAnaElemTool>& anatools)
    : G4UserTrackingAction(),
      m_anaelemtools(anatools) {

}

TrackingAction::~TrackingAction() {

}

void
TrackingAction::PreUserTrackingAction(const G4Track* aTrack) {

}

void
TrackingAction::PostUserTrackingAction(const G4Track* aTrack) {

}


