#ifndef EventAction_h
#define EventAction_h

#include <GaudiKernel/ToolHandle.h>

#include <DetSimInterface/IAnaElemTool.h>

#include "G4UserEventAction.hh"

class G4Event;

class EventAction: public G4UserEventAction {
public:

    EventAction(ToolHandleArray<IAnaElemTool>&);
    ~EventAction();

    void BeginOfEventAction(const G4Event*) override;
    void EndOfEventAction(const G4Event*) override;

private:
    ToolHandleArray<IAnaElemTool>& m_anaelemtools;
};

#endif
