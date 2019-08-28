#include "HepevtRdr.h"
#include "GenEvent.h"

#include "lcio.h"  //LCIO
#include "EVENT/LCIO.h"
#include "LCAscHepRdr.h"
#include "IMPL/MCParticleImpl.h"


#include "plcio/MCParticle.h" //plcio
#include "plcio/MCParticleObj.h"
#include "plcio/MCParticleCollection.h"
#include "plcio/DoubleThree.h"
#include "plcio/FloatThree.h"
#include "plcio/EventHeaderCollection.h"



#include <iostream>
#include <vector>
#include <fstream>


using namespace lcio;
using namespace IMPL;
using namespace plcio;
using namespace std;

typedef enum HEPFILEFORMATS
{
  stdhep = 0,
  HEPEvt,
  hepevt,
  slcio
}  HEPFILEFORMAT;


HepevtRdr::HepevtRdr(string name){

m_hepevt_rdr = new UTIL::LCAscHepRdr(name.c_str(), hepevt);
m_processed_event=0;
std::cout<<"initial hepevt_rdr"<<std::endl;
}

HepevtRdr::~HepevtRdr(){
delete m_hepevt_rdr;
}

bool HepevtRdr::mutate(MyHepMC::GenEvent& event){
    LCCollectionVec* mc_vec = m_hepevt_rdr->readEvent();
    if(mc_vec==nullptr) return false;
    m_processed_event ++;
    int n_mc = mc_vec->getNumberOfElements();
    std::cout<<"Read event :"<< m_processed_event <<", mc size :"<< n_mc <<std::endl;
    std::map<int, int> pmcid_lmcid;
    for (int i=0; i < n_mc; i++){
        MCParticleImpl* mc = (MCParticleImpl*) mc_vec->getElementAt(i);
        //std::cout<<"At mc :"<< i <<std::endl;
        plcio::MCParticle mcp = event.m_mc_vec.create();
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
        mcp.setMomentum           (FloatThree(float(mc->getMomentum()[0]), float(mc->getMomentum()[1]), float(mc->getMomentum()[2]) ));
        mcp.setMomentumAtEndpoint (FloatThree(float(mc->getMomentumAtEndpoint()[0]), float(mc->getMomentumAtEndpoint()[1]), float(mc->getMomentumAtEndpoint()[2]) ));
        mcp.setSpin               (mc->getSpin());
        mcp.setColorFlow          (mc->getColorFlow());
    }
    // second loop for setting parents and daughters
    
    for (int i=0; i < n_mc; i++){
        MCParticleImpl* mc = (MCParticleImpl*) mc_vec->getElementAt(i);
        const MCParticleVec & mc_parents = mc->getParents();
        const MCParticleVec & mc_daughters = mc->getDaughters();
        plcio::MCParticle pmc = event.m_mc_vec.at(i);
        //std::cout<<"mc at "<< i<<", parent size "<<mc_parents.size() <<std::endl;
        for(unsigned int j=0; j< mc_parents.size(); j++){int p_id = mc_parents.at(j)->id();
                                                 //std::cout<<"parent id "<<p_id<<std::endl;
                                                 pmc.addParent( event.m_mc_vec.at( pmcid_lmcid.at(p_id) ) );
                                                }
        //std::cout<<"mc at "<< i<<", daughter size "<<mc_daughters.size() <<std::endl;
        for(unsigned int j=0; j< mc_daughters.size(); j++){int d_id = mc_daughters.at(j)->id();
                                                 //std::cout<<"daughter id "<<d_id<<std::endl;
                                                 pmc.addDaughter( event.m_mc_vec.at( pmcid_lmcid.at(d_id) ) );
                                                }
    }
     
    event.SetEventHeader( m_processed_event, -99, 9999, "Generator");
    //std::cout<<"end event :"<< m_processed_event <<std::endl;
    return true;
}

bool HepevtRdr::isEnd(){
return false;
}

bool HepevtRdr::configure(){
return true;
}

bool HepevtRdr::finish(){
return true;
}
