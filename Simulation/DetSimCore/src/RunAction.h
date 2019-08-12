#ifndef RunAction_h
#define RunAction_h

#include "G4UserRunAction.hh"

class G4Run;


class RunAction: public G4UserRunAction {

public:
    RunAction();
    ~RunAction();

    virtual void BeginOfRunAction(const G4Run*);
    virtual void EndOfRunAction(const G4Run*);

};

#endif
