#ifndef IDetElemTool_h
#define IDetElemTool_h

// IDetElemTool is used to wrap the construction of G4LogicalVolume.
// Please note that the placement of logical volume is fixed in the code.
// If necessary, another IDetElemPosTool can be used to produce the positions
// of the daughters.
// IDetElemTool should represent the high level detectors/modules.

#include "GaudiKernel/IAlgTool.h"

class G4LogicalVolume;

class IDetElemTool: virtual public IAlgTool {
public:
    DeclareInterfaceID(IDetElemTool, 0, 1);

    virtual ~IDetElemTool() {}

    // return the constructed detector
    virtual G4LogicalVolume* getLV() = 0;
    virtual void ConstructSDandField() {}
};

#endif
