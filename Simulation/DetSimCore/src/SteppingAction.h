#ifndef SteppingAction_h
#define SteppingAction_h

#include <GaudiKernel/ToolHandle.h>

#include <DetSimInterface/IAnaElemTool.h>

#include "G4UserSteppingAction.hh"

class G4Step;

class SteppingAction: public G4UserSteppingAction {

public:
    SteppingAction(ToolHandleArray<IAnaElemTool>&);
    ~SteppingAction();

    void UserSteppingAction(const G4Step*) override;

private:
    ToolHandleArray<IAnaElemTool>& m_anaelemtools;
};

#endif
