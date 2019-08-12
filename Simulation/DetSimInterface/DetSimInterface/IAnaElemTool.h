#ifndef IAnaElemTool_h
#define IAnaElemTool_h

#include "GaudiKernel/IAlgTool.h"

class G4Run;
class G4Event;
class G4Track;
class G4Step;

#include "G4ClassificationOfNewTrack.hh"

class IAnaElemTool : virtual public IAlgTool {
public:
    DeclareInterfaceID(IAnaElemTool, 0, 1);

    virtual ~IAnaElemTool() {}

    // Run
    virtual void BeginOfRunAction(const G4Run*) {}
    virtual void EndOfRunAction(const G4Run*) {}

    // Event
    virtual void BeginOfEventAction(const G4Event*) {}
    virtual void EndOfEventAction(const G4Event*) {}

    // Stacking
    virtual G4ClassificationOfNewTrack ClassifyNewTrack(const G4Track*) {return fUrgent;}
    virtual void NewStage() {}
    virtual void PrepareNewEvent() {}

    // Tracking
    virtual void PreUserTrackingAction(const G4Track*) {}
    virtual void PostUserTrackingAction(const G4Track*) {}

    // Stepping
    virtual void UserSteppingAction(const G4Step*) {}

};


#endif
