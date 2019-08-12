#ifndef DetSimAlg_h
#define DetSimAlg_h

#include <string>
#include <vector>

#include <GaudiKernel/Algorithm.h>
#include <GaudiKernel/Property.h>
#include <GaudiKernel/ToolHandle.h>

#include <DetSimInterface/IDetSimSvc.h>
#include <DetSimInterface/IAnaElemTool.h>

class DetSimAlg: public Algorithm {
public:
    DetSimAlg(const std::string& name, ISvcLocator* pSvcLocator);

    StatusCode initialize() override;
    StatusCode execute() override;
    StatusCode finalize() override;

private:
    SmartIF<IDetSimSvc> m_detsimsvc;
    ToolHandleArray<IAnaElemTool> m_anaelemtools;

private:

    Gaudi::Property<std::vector<std::string>> m_run_macs{this, "RunMacs"};
    Gaudi::Property<std::vector<std::string>> m_run_cmds{this, "RunCmds"};
    Gaudi::Property<std::vector<std::string>> m_vis_macs{this, "VisMacs"};

    Gaudi::Property<std::string> m_physics_lists_name{this, "PhysicsList", "QGSP_BERT"};

    Gaudi::Property<std::vector<std::string>> m_ana_elems{this, "AnaElems"};


private:
    int i_event;
};



#endif
