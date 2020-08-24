#include "StdHepRdr.h"
#include "GenEvent.h"

#include "lcio.h"  //LCIO
#include "EVENT/LCIO.h"
#include "UTIL/LCStdHepRdrNew.h"
#include "IMPL/MCParticleImpl.h"


#include "edm4hep/MCParticle.h" //edm4hep
#include "edm4hep/MCParticleObj.h"
#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/EventHeaderCollection.h"



#include <iostream>
#include <vector>
#include <fstream>


using namespace lcio;
using namespace IMPL;
using namespace edm4hep;
using namespace std;

DECLARE_COMPONENT(StdHepRdr)

StdHepRdr::~StdHepRdr(){
    delete m_stdhep_rdr;
}

bool StdHepRdr::mutate(MyHepMC::GenEvent& event){
    if(isEnd()) return false;
    LCCollectionVec* mc_vec = m_stdhep_rdr->readEvent();
    m_processed_event ++;
    int n_mc = mc_vec->getNumberOfElements();
    //std::cout<<"Debug: Read event :"<< m_processed_event <<", mc size :"<< n_mc <<std::endl;
    std::map<int, int> pmcid_lmcid;
    for (int i=0; i < n_mc; i++){
        MCParticleImpl* mc = (MCParticleImpl*) mc_vec->getElementAt(i);
        //std::cout<<"At mc :"<< i <<std::endl;
        edm4hep::MCParticle mcp = event.m_mc_vec.create();
        pmcid_lmcid.insert(std::pair<int, int>(mc->id(),i));
        //std::cout<<"map<id,i>:"<<mc->id()<<","<< i <<std::endl;
                                 
        mcp.setPDG                (mc->getPDG());  
        mcp.setGeneratorStatus    (mc->getGeneratorStatus());
        mcp.setSimulatorStatus    (mc->getSimulatorStatus());
        mcp.setCharge             (mc->getCharge());
        mcp.setTime               (mc->getTime());
        mcp.setMass               (mc->getMass());
        mcp.setVertex             (mc->getVertex()); 
        mcp.setEndpoint           (mc->getEndpoint());
        mcp.setMomentum           (Vector3f(float(mc->getMomentum()[0]), float(mc->getMomentum()[1]), float(mc->getMomentum()[2]) ));
        mcp.setMomentumAtEndpoint (Vector3f(float(mc->getMomentumAtEndpoint()[0]), float(mc->getMomentumAtEndpoint()[1]), float(mc->getMomentumAtEndpoint()[2]) ));
        mcp.setSpin               (mc->getSpin());
        mcp.setColorFlow          (mc->getColorFlow());
    }
    // second loop for setting parents and daughters
    
    for (int i=0; i < n_mc; i++){
        MCParticleImpl* mc = (MCParticleImpl*) mc_vec->getElementAt(i);
        const MCParticleVec & mc_parents = mc->getParents();
        const MCParticleVec & mc_daughters = mc->getDaughters();
        edm4hep::MCParticle pmc = event.m_mc_vec.at(i);
        //std::cout<<"mc at "<< i<<", parent size "<<mc_parents.size() <<std::endl;
        for(unsigned int j=0; j< mc_parents.size(); j++){int p_id = mc_parents.at(j)->id();
                                                 //std::cout<<"parent id "<<p_id<<std::endl;
                                                 pmc.addToParents( event.m_mc_vec.at( pmcid_lmcid.at(p_id) ) );
                                                }
        //std::cout<<"mc at "<< i<<", daughter size "<<mc_daughters.size() <<std::endl;
        for(unsigned int j=0; j< mc_daughters.size(); j++){int d_id = mc_daughters.at(j)->id();
                                                 //std::cout<<"daughter id "<<d_id<<std::endl;
                                                 pmc.addToDaughters( event.m_mc_vec.at( pmcid_lmcid.at(d_id) ) );
                                                }
    }
     
    event.SetEventHeader( m_processed_event, -99, 9999, "Generator");
    //std::cout<<"end event :"<< m_processed_event <<std::endl;

    delete mc_vec;

    return true;
}

bool StdHepRdr::isEnd(){
if(m_processed_event == m_total_event) {std::cout<<"Have read all events, end now."<<std::endl; return true;}
else return false;
}

bool StdHepRdr::configure_gentool(){
    m_stdhep_rdr = new LCStdHepRdrNew(m_filename.value().c_str());
    m_stdhep_rdr->printHeader();
    if (m_stdhep_rdr->getNumberOfEvents()<=1) {
        return false;
    }

    m_total_event = m_stdhep_rdr->getNumberOfEvents() - 1 ;
    m_processed_event=0;

    return true;
}

bool StdHepRdr::finish(){
    return true;
}

StatusCode
StdHepRdr::initialize() {
    StatusCode sc;
    if (not configure_gentool()) {
        error() << "failed to initialize." << endmsg;
        return StatusCode::FAILURE;
    }

    return sc;
}

StatusCode
StdHepRdr::finalize() {
    StatusCode sc;
    if (not finish()) {
        error() << "Failed to finalize." << endmsg;
        return StatusCode::FAILURE;
    }
    return sc;
}
