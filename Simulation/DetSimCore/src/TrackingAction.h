#ifndef TrackingAction_h
#define TrackingAction_h

#include <GaudiKernel/ToolHandle.h>

#include <DetSimInterface/IAnaElemTool.h>

#include "G4UserTrackingAction.hh"

class TrackingAction: public G4UserTrackingAction {

public:

    TrackingAction(ToolHandleArray<IAnaElemTool>&);
    ~TrackingAction();

    void PreUserTrackingAction(const G4Track*);
    void PostUserTrackingAction(const G4Track*);

private:
    ToolHandleArray<IAnaElemTool>& m_anaelemtools;
    
};

#endif
