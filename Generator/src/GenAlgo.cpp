#include "GenAlgo.h"

#include "GaudiKernel/IEventProcessor.h"
#include "GaudiKernel/IAppMgrUI.h"
#include "GaudiKernel/GaudiException.h"


#include "edm4hep/MCParticleCollection.h"//plico

#include <iostream>
#include <vector>
#include <fstream>

#include "IGenTool.h"
#include "GenEvent.h"


DECLARE_COMPONENT(GenAlgo)

GenAlgo::GenAlgo(const std::string& name, ISvcLocator* pSvcLocator): GaudiAlgorithm(name, pSvcLocator) {
    declareProperty("MCParticleGen", m_hdl, "MCParticle collection (at Generator phase)");
    declareProperty("GenTools", m_genToolNames, "List of GenTools");
    m_evtid = 0;

}

StatusCode
GenAlgo::initialize() {
    if (m_genToolNames.size()==0) {
        error() << "Please specify the gentools to be used in GenAlgo." << endmsg;
        return StatusCode::FAILURE;
    }

    for (auto gtname: m_genToolNames) {
        m_genTools.push_back(gtname);
    }
    
    return StatusCode::SUCCESS;

}

StatusCode
GenAlgo::execute() {
    m_evtid++;
    auto mcCol = m_hdl.createAndPut();
    MyHepMC::GenEvent m_event(*mcCol);

    for(auto gentool: m_genTools) {
        if (gentool->mutate(m_event)) {} 
        else {
            warning() << "Have read all events, stop now." << endmsg; 
            auto ep = serviceLocator()->as<IEventProcessor>();
            if ( !ep ) {
                error() << "Cannot get IEventProcessor" << endmsg;
                return StatusCode::FAILURE;
            }
            ep->stopRun();
            return StatusCode::SUCCESS;
            
        }
    }

    return StatusCode::SUCCESS;

}

StatusCode
GenAlgo::finalize() {
    return StatusCode::SUCCESS;
}
