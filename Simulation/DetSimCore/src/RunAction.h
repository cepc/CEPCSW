#ifndef RunAction_h
#define RunAction_h

#include <GaudiKernel/ToolHandle.h>

#include <DetSimInterface/IAnaElemTool.h>

#include "G4UserRunAction.hh"

class G4Run;


class RunAction: public G4UserRunAction {

public:
    RunAction(ToolHandleArray<IAnaElemTool>&);
    ~RunAction();

    virtual void BeginOfRunAction(const G4Run*);
    virtual void EndOfRunAction(const G4Run*);
private:
    ToolHandleArray<IAnaElemTool>& m_anaelemtools;

};

#endif
