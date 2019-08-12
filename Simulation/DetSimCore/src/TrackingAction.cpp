#include "TrackingAction.h"

TrackingAction::TrackingAction(ToolHandleArray<IAnaElemTool>& anatools)
    : G4UserTrackingAction(),
      m_anaelemtools(anatools) {

}

TrackingAction::~TrackingAction() {

}

void
TrackingAction::PreUserTrackingAction(const G4Track* aTrack) {
    for (auto ana: m_anaelemtools) {
        ana->PreUserTrackingAction(aTrack);
    }
}

void
TrackingAction::PostUserTrackingAction(const G4Track* aTrack) {
    for (auto ana: m_anaelemtools) {
        ana->PostUserTrackingAction(aTrack);
    }
}


