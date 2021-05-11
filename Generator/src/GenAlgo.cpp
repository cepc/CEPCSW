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
// #include "StdHepRdr.h"
// #include "HepevtRdr.h"// not correct still
// #include "SLCIORdr.h"
// #include "HepMCRdr.h"
// #include "GenPrinter.h"
// #include "GenWriter.h"

using namespace std;

DECLARE_COMPONENT(GenAlgo)

GenAlgo::GenAlgo(const std::string& name, ISvcLocator* pSvcLocator): GaudiAlgorithm(name, pSvcLocator) {
    declareProperty("MCParticle", m_hdl, "MCParticle collection (output)");
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
    
    // cout << "initialize start" << endl; 
    // string generatorName = m_input_file.value();
    // string outputName    = m_output_file.value();
    // string format        = m_input_format.value();
    // IGenTool* gen_reader;
    // if(format=="stdhep") gen_reader  = new StdHepRdr(generatorName);    
    // else if(format=="slcio") gen_reader  = new SLCIORdr(generatorName);    
    // else if(format=="hepmc") gen_reader  = new HepMCRdr(generatorName);    
    // else{cout << "Error : unsupport format for generator input file" << endl; return StatusCode::FAILURE; }
    // //IGenTool* gen_reader  = new HepevtRdr(generatorName);    
    // m_genTools.push_back(gen_reader);
    // if(m_print.value()) {
    //     IGenTool* gen_printer = new GenPrinter(generatorName);    
    //     m_genTools.push_back(gen_printer);
    // }
    //IGenTool* gen_writer  = new GenWriter (outputName);    
    //m_genTools.push_back(gen_writer);

    // cout << "initialize done" << endl; 
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
            cout << "Have read all events, stop now." << endl; 
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
    // cout << "finalize" << endl; 
    // for(auto gentool: m_genTools) {
    //     if (gentool->finish()) {} 
    //     else {cout << "finish Failed" << endl; return StatusCode::FAILURE; }
    // }
    return StatusCode::SUCCESS;
}
