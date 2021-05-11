#include "SLCIORdr.h"
#include "GenEvent.h"

#include "lcio.h"  //LCIO
#include "LCIOSTLTypes.h"
#include "IOIMPL/LCFactory.h"
#include "EVENT/LCIO.h"
#include "EVENT/LCEvent.h"
#include "EVENT/LCCollection.h"
#include "IO/LCReader.h"
#include "IMPL/MCParticleImpl.h"
#include "IMPL/LCCollectionVec.h"


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

DECLARE_COMPONENT(SLCIORdr)

SLCIORdr::~SLCIORdr(){
delete m_slcio_rdr;
}

bool SLCIORdr::mutate(MyHepMC::GenEvent& event){

      IMPL::LCCollectionVec* lcMCVec = new IMPL::LCCollectionVec(LCIO::MCPARTICLE);
      EVENT::LCEvent *lcEvent = m_slcio_rdr->readNextEvent(LCIO::UPDATE);
      LCCollection *lcCol = NULL;
      if(lcEvent) lcCol = lcEvent->getCollection("MCParticle");
      else return false;
      if(lcCol){
	MCParticleImpl* p;
	//MCParticleImpl* d;
	int NHEP = lcCol->getNumberOfElements();
	lcMCVec->parameters() = lcCol->parameters();
	lcMCVec->resize(NHEP);
	for( int IHEP=0; IHEP<NHEP; IHEP++ ){
	  MCParticleImpl* mcp = new MCParticleImpl();
	  EVENT::MCParticle* in = dynamic_cast<EVENT::MCParticle*>(lcCol->getElementAt(IHEP));
	  lcMCVec->at(IHEP)  = mcp ;
	  mcp->setPDG(in->getPDG());
	  mcp->setCharge(in->getCharge()) ;
	  mcp->setMomentum(in->getMomentum());
	  mcp->setMass(in->getMass());
	  mcp->setVertex(in->getVertex());
	  mcp->setGeneratorStatus(in->getGeneratorStatus());
	  mcp->setSimulatorStatus(in->getSimulatorStatus());
	  mcp->setTime(in->getTime());
	  mcp->setSpin(in->getSpin());
	}
	for( int IHEP=0; IHEP<NHEP; IHEP++ ){
	  lcio::MCParticleImpl* mcp = dynamic_cast<MCParticleImpl*>(lcMCVec->getElementAt(IHEP));
	  EVENT::MCParticleVec parents = (dynamic_cast<EVENT::MCParticle*>(lcCol->getElementAt(IHEP)))->getParents();
	  int np = parents.size();
	  //cout << "DEBUG: the " << IHEP << "th particle has " << np << " parents ";
	  int i;
	  for(int ip=0;ip<np;ip++){
	    p = dynamic_cast<MCParticleImpl*>(parents[ip]);
	    for(i=0;i<NHEP;i++){
	      if(p==lcCol->getElementAt(i)) break;
	    }
	    if(i==NHEP) cout << "Heedm4hepInterfaceNew: error" << endl;
	    mcp->addParent(dynamic_cast<MCParticleImpl*>(lcMCVec->getElementAt(i)));
	  }
	  //ignore daughter table, auto-regonize relationship while addParent()
          /*
	  EVENT::MCParticleVec daughters = (dynamic_cast<EVENT::MCParticle*>(lcCol->getElementAt(IHEP)))->getDaughters();
          int nd = daughters.size();
	  //cout << "and " << nd << " daughters." << endl;
          for(int id=0;id<nd;id++){
            d = dynamic_cast<MCParticleImpl*>(daughters[id]);
            for(i=0;i<NHEP;i++){
              if(d==lcCol->getElementAt(i)) break;
            }
            if(i==NHEP) cout << "Heedm4hepInterfaceNew: error" << endl;
            mcp->addToDaughters(dynamic_cast<MCParticleImpl*>(lcMCVec->getElementAt(i)));
          }
          */
          
	}
      }
      else{
	cout << "Debug: no MCParticle Collection is read!" << endl;
        return false;
      }


    m_processed_event ++;
    int n_mc = lcMCVec->getNumberOfElements();
    std::cout<<"Read event :"<< m_processed_event <<", mc size :"<< n_mc <<std::endl;
    std::map<int, int> pmcid_lmcid;
    for (int i=0; i < n_mc; i++){
        MCParticleImpl* mc = (MCParticleImpl*) lcMCVec->getElementAt(i);
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
        MCParticleImpl* mc = (MCParticleImpl*) lcMCVec->getElementAt(i);
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
    //for(unsigned int i=0; i< lcMCVec->size(); i++) delete lcMCVec->at(i);
    delete lcMCVec;  
    return true;
}

bool SLCIORdr::isEnd(){
return false;
}

bool SLCIORdr::configure_gentool(){
    m_slcio_rdr = IOIMPL::LCFactory::getInstance()->createLCReader();
    m_slcio_rdr->open(m_filename.value().c_str());
    m_processed_event=0;


    return true;
}

bool SLCIORdr::finish(){
return true;
}


StatusCode
SLCIORdr::initialize() {
    StatusCode sc;
    if (not configure_gentool()) {
        error() << "failed to initialize." << endmsg;
        return StatusCode::FAILURE;
    }

    return sc;
}

StatusCode
SLCIORdr::finalize() {
    StatusCode sc;
    if (not finish()) {
        error() << "Failed to finalize." << endmsg;
        return StatusCode::FAILURE;
    }
    return sc;
}
