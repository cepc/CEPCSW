#ifndef WorldDetElemTool_h
#define WorldDetElemTool_h

#include "GaudiKernel/AlgTool.h"
#include "DetSimInterface/IDetElemTool.h"

class WorldDetElemTool: public extends<AlgTool, IDetElemTool> {

public:
    using extends::extends;

    G4LogicalVolume* getLV() override;
    void ConstructSDandField() override;

    StatusCode initialize() override;
    StatusCode finalize() override;

private:
    double m_x;
    double m_y;
    double m_z;
};

#endif
