#ifndef GenAlgo_h
#define GenAlgo_h


#include <GaudiKernel/Algorithm.h>
#include <Gaudi/Property.h>
#include <GaudiKernel/ToolHandle.h>

#include "GaudiAlg/GaudiAlgorithm.h"
#include "k4FWCore/DataHandle.h"

#include "GenEvent.h"

class IGenTool;
namespace plcio {
    class MCParticleCollection;
}


using namespace std;

class GenAlgo: public GaudiAlgorithm {

public:
    GenAlgo(const std::string& name, ISvcLocator* pSvcLocator);

    StatusCode initialize() override;
    StatusCode execute() override;
    StatusCode finalize() override;

private:
    Gaudi::Property<std::string> m_input_file{this, "Input", "NULL"};
    Gaudi::Property<std::string> m_input_format{this, "FileFormat", "NULL"};
    Gaudi::Property<std::string> m_output_file{this, "OutputRootFile", "NULL"};
    Gaudi::Property<bool> m_print{this, "PrintEvent", "NULL"};
    Gaudi::Property<bool> m_do_write{this, "WriteFile", "NULL"};

    std::vector<std::string> m_genToolNames;                                                         
    // std::vector<IGenTool*> m_genTools;
    ToolHandleArray<IGenTool> m_genTools;

    int m_evtid;                               
    int m_evtMax;
    //MyHepMC::GenEvent m_event;
    DataHandle<edm4hep::MCParticleCollection> m_hdl{"MCParticle", Gaudi::DataHandle::Writer, this};


};


#endif
