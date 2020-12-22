#ifndef AnExampleDetElemTool_h
#define AnExampleDetElemTool_h

#include "GaudiKernel/AlgTool.h"
#include <Gaudi/Property.h>
#include <GaudiKernel/ToolHandle.h>

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

#include "DetInterface/IGeomSvc.h"
#include "DetSimInterface/IDetElemTool.h"
#include "DetSimInterface/ISensDetTool.h"


class AnExampleDetElemTool: public extends<AlgTool, IDetElemTool> {

public:
    using extends::extends;

    G4LogicalVolume* getLV() override;
    void ConstructSDandField() override;

    StatusCode initialize() override;
    StatusCode finalize() override;

private:
    Gaudi::Property<double> m_x{this, "X", 30.*m};
    Gaudi::Property<double> m_y{this, "Y", 30.*m};
    Gaudi::Property<double> m_z{this, "Z", 30.*m};
    // DD4hep XML compact file path
    Gaudi::Property<std::string> m_dd4hep_xmls{this, "detxml"};

    SmartIF<IGeomSvc> m_geosvc;
    ToolHandle<ISensDetTool> m_calo_sdtool;
    ToolHandle<ISensDetTool> m_driftchamber_sdtool;
    ToolHandle<ISensDetTool> m_tpc_sdtool;
};

#endif
