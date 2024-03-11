#include "HepevtRdr.h"
#include "GenEvent.h"

#include "lcio.h"  //LCIO
#include "EVENT/LCIO.h"
#include "LCAscHepRdr.h"
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

typedef enum HEPFILEFORMATS
{
  stdhep = 0,
  HEPEvt,
  hepevt,
  slcio
}  HEPFILEFORMAT;


DECLARE_COMPONENT(HepevtRdr)

HepevtRdr::~HepevtRdr(){
    delete m_hepevt_rdr;
}

StatusCode HepevtRdr::initialize() {
    StatusCode sc;
    if (not configure_gentool()) {
        error() << "failed to initialize." << endmsg;
        return StatusCode::FAILURE;
    }

    return sc;
}

StatusCode HepevtRdr::finalize() {
    StatusCode sc;
    if (not finish()) {
        error() << "Failed to finalize." << endmsg;
        return StatusCode::FAILURE;
    }
    return sc;
}


bool HepevtRdr::configure_gentool(){
    int format = hepevt;
    if (m_format == "HEPEvt") {
        format = HEPEvt;
    } else if (m_format == "hepevt") {
        format = hepevt;
    }

    m_hepevt_rdr = new UTIL::LCAscHepRdr(m_filename.value().c_str(), format);
    m_processed_event=0;
    std::cout<<"initial hepevt_rdr"<<std::endl;
    return true;
}

bool HepevtRdr::mutate(MyHepMC::GenEvent& event){
    LCCollectionVec* mc_vec = m_hepevt_rdr->readEvent();
    if(mc_vec==nullptr) return false;
    m_processed_event ++;
    int n_mc = mc_vec->size();
    std::cout<<"Read event :"<< m_processed_event <<", mc size :"<< n_mc <<std::endl;
    std::map<int, int> pmcid_lmcid;
    for (int i=0; i < n_mc; i++){
        MCParticleImpl* mc = (MCParticleImpl*) mc_vec->getElementAt(i);
        // std::cout<<"At mc :"<< i <<std::endl;
        auto mcp = event.m_mc_vec.create();
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
        mcp.setMomentum           (edm4hep::Vector3f(float(mc->getMomentum()[0]), float(mc->getMomentum()[1]), float(mc->getMomentum()[2]) ));
        mcp.setMomentumAtEndpoint (edm4hep::Vector3f(float(mc->getMomentumAtEndpoint()[0]), float(mc->getMomentumAtEndpoint()[1]), float(mc->getMomentumAtEndpoint()[2]) ));
        mcp.setSpin               (mc->getSpin());
        mcp.setColorFlow          (mc->getColorFlow());
    }
    // second loop for setting parents and daughters
    
    for (int i=0; i < n_mc; i++){
        MCParticleImpl* mc = (MCParticleImpl*) mc_vec->getElementAt(i);
        const MCParticleVec & mc_parents = mc->getParents();
        const MCParticleVec & mc_daughters = mc->getDaughters();
        auto pmc = event.m_mc_vec.at(i);
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
    return true;
}

bool HepevtRdr::isEnd(){
    return false;
}

bool HepevtRdr::finish(){
    return true;
}
