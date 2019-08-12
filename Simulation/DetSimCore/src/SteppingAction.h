#ifndef SteppingAction_h
#define SteppingAction_h

#include "G4UserSteppingAction.hh"

class G4Step;

class SteppingAction: public G4UserSteppingAction {

public:
    SteppingAction();
    ~SteppingAction();

    void UserSteppingAction(const G4Step*) override;

};

#endif
