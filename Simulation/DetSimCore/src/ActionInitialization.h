#ifndef ActionInitialization_h
#define ActionInitialization_h

#include <GaudiKernel/ToolHandle.h>

#include <DetSimInterface/IAnaElemTool.h>


#include "G4VUserActionInitialization.hh"

class ActionInitialization: public G4VUserActionInitialization {
public:

    ActionInitialization(ToolHandleArray<IAnaElemTool>&);
    ~ActionInitialization();

    void BuildForMaster() const override;
    void Build() const override;

private:
    ToolHandleArray<IAnaElemTool>& m_anaelemtools;
};


#endif
