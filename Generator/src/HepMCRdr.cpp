#include "HepMCRdr.h"
#include "GenEvent.h"

#include "HepMC/IO_GenEvent.h"//HepMC
#include "HepMC/GenEvent.h"
#include "HepMC/GenVertex.h"
#include "HepMC/Polarization.h"


#include "plcio/MCParticle.h" //plcio
#include "plcio/MCParticleObj.h"
#include "plcio/MCParticleCollection.h"
#include "plcio/DoubleThree.h"
#include "plcio/FloatThree.h"
#include "plcio/EventHeaderCollection.h"



#include <iostream>
#include <vector>
#include <fstream>


using namespace plcio;
using namespace std;


HepMCRdr::HepMCRdr(string name){

ascii_in = new HepMC::IO_GenEvent(name.c_str(),std::ios::in);

m_processed_event=0;
}

HepMCRdr::~HepMCRdr(){
delete ascii_in;
}

bool HepMCRdr::mutate(MyHepMC::GenEvent& event){

    HepMC::GenEvent* evt = ascii_in->read_next_event();
    if(!evt) return false;
    m_processed_event ++;
    int n_mc = evt->particles_size();
    //std::cout<<"Read event :"<< m_processed_event <<", mc size :"<< n_mc <<std::endl;
    std::map<int, int> pmcid_lmcid;
    int index = 0 ;
    for ( HepMC::GenEvent::particle_iterator p = evt->particles_begin(); p != evt->particles_end(); ++p ) {
        //std::cout<<"start mc "<<index<<std::endl;
        plcio::MCParticle mcp = event.m_mc_vec.create();
        pmcid_lmcid.insert(std::pair<int, int>((*p)->barcode(),index));
        index++;
        //std::cout<<"map<id,i>:"<<mc->id()<<","<< i <<std::endl;
                                 
        mcp.setPDG                ((*p)->pdg_id());  
        mcp.setGeneratorStatus    ((*p)->status());
        mcp.setSimulatorStatus    (999);
        mcp.setCharge             (999);
        mcp.setTime               (999);
        mcp.setMass               ((*p)->generated_mass());
        if ( (*p)->production_vertex() ){
            HepMC::GenVertex* vertex_pro =  (*p)->production_vertex();
            double three[3] = {vertex_pro->point3d().x(), vertex_pro->point3d().y(), vertex_pro->point3d().z()};
            mcp.setVertex             (plcio::DoubleThree (three)); 
        }
        else mcp.setVertex(plcio::DoubleThree()); 
        if ( (*p)->end_vertex() ){
            HepMC::GenVertex* vertex_end =  (*p)->end_vertex();
            double three[3] = {vertex_end->point3d().x(), vertex_end->point3d().y(), vertex_end->point3d().z()};
            mcp.setEndpoint           (plcio::DoubleThree (three));
        } 
        else mcp.setEndpoint (plcio::DoubleThree());
        mcp.setMomentum           (plcio::FloatThree(float((*p)->momentum().px()), float((*p)->momentum().py()), float((*p)->momentum().pz()) ));
        mcp.setMomentumAtEndpoint (plcio::FloatThree(float((*p)->momentum().px()), float((*p)->momentum().py()), float((*p)->momentum().pz()) ));
        const HepMC::Polarization & polar = (*p)->polarization();
        mcp.setSpin               (plcio::FloatThree(polar.normal3d().x(), polar.normal3d().y(), polar.normal3d().z()) );
        int two[2] = {1, (*p)->flow(1)};
        mcp.setColorFlow          (plcio::IntTwo (two) );
    }
    // second loop for setting parents and daughters
    index = 0 ;
    for ( HepMC::GenEvent::particle_iterator p = evt->particles_begin(); p != evt->particles_end(); ++p ) {
        plcio::MCParticle pmc = event.m_mc_vec.at(index);
        index++;
        if ( (*p)->production_vertex() ) {
            for ( HepMC::GenVertex::particle_iterator mother = (*p)->production_vertex()-> particles_begin(HepMC::parents); mother != (*p)->production_vertex()-> particles_end(HepMC::parents); ++mother ) {
                pmc.addParent( event.m_mc_vec.at( pmcid_lmcid.at((*mother)->barcode()) ) );
            }
        }
        if ( (*p)->end_vertex() ) {
            for ( HepMC::GenVertex::particle_iterator des =(*p)->end_vertex()-> particles_begin(HepMC::descendants); des != (*p)->end_vertex()-> particles_end(HepMC::descendants); ++des ) {
                pmc.addDaughter( event.m_mc_vec.at( pmcid_lmcid.at((*des)->barcode()) ) );
                }
        }   
    }
    
    event.SetEventHeader( m_processed_event, -99, 9999, "Generator");
    //std::cout<<"end event :"<< m_processed_event <<std::endl;
    delete evt;
    return true;
}

bool HepMCRdr::isEnd(){
return false;
}

bool HepMCRdr::configure(){
return true;
}

bool HepMCRdr::finish(){
return true;
}
