#ifndef TrackingAction_h
#define TrackingAction_h

#include "G4UserTrackingAction.hh"

class TrackingAction: public G4UserTrackingAction {

public:

    TrackingAction();
    ~TrackingAction();

    void PreUserTrackingAction(const G4Track*);
    void PostUserTrackingAction(const G4Track*);



};

#endif
