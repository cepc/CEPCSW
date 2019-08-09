#ifndef DetSimAlg_h
#define DetSimAlg_h

#include <string>
#include <vector>

#include <GaudiKernel/Algorithm.h>
#include <GaudiKernel/Property.h>

#include <DetSimInterface/IDetSimSvc.h>

class DetSimAlg: public Algorithm {
public:
    DetSimAlg(const std::string& name, ISvcLocator* pSvcLocator);

    StatusCode initialize() override;
    StatusCode execute() override;
    StatusCode finalize() override;

private:
    SmartIF<IDetSimSvc> m_detsimsvc;

private:

    Gaudi::Property<std::vector<std::string>> m_run_macs{this, "RunMacs"};
    Gaudi::Property<std::vector<std::string>> m_run_cmds{this, "RunCmds"};
    Gaudi::Property<std::vector<std::string>> m_vis_macs{this, "VisMacs"};

    Gaudi::Property<std::string> m_physics_lists_name{this, "PhysicsList", "QGSP_BERT"};
private:
    int i_event;
};



#endif
