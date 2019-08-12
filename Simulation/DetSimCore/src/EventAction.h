#ifndef EventAction_h
#define EventAction_h

#include "G4UserEventAction.hh"

class G4Event;

class EventAction: public G4UserEventAction {
public:

    EventAction();
    ~EventAction();

    void BeginOfEventAction(const G4Event*) override;
    void EndOfEventAction(const G4Event*) override;

};

#endif
