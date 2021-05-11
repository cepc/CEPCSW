#include "HepMCRdr.h"
#include "GenEvent.h"

#include "HepMC/IO_GenEvent.h"//HepMC
#include "HepMC/GenEvent.h"
#include "HepMC/GenVertex.h"
#include "HepMC/Polarization.h"


#include "edm4hep/MCParticle.h" //edm4hep
#include "edm4hep/MCParticleObj.h"
#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/EventHeaderCollection.h"



#include <iostream>
#include <vector>
#include <fstream>


using namespace edm4hep;
using namespace std;

DECLARE_COMPONENT(HepMCRdr)

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
        edm4hep::MCParticle mcp = event.m_mc_vec.create();
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
            mcp.setVertex             (edm4hep::Vector3d (three)); 
        }
        else mcp.setVertex(edm4hep::Vector3d()); 
        if ( (*p)->end_vertex() ){
            HepMC::GenVertex* vertex_end =  (*p)->end_vertex();
            double three[3] = {vertex_end->point3d().x(), vertex_end->point3d().y(), vertex_end->point3d().z()};
            mcp.setEndpoint           (edm4hep::Vector3d (three));
        } 
        else mcp.setEndpoint (edm4hep::Vector3d());
        mcp.setMomentum           (edm4hep::Vector3f(float((*p)->momentum().px()), float((*p)->momentum().py()), float((*p)->momentum().pz()) ));
        mcp.setMomentumAtEndpoint (edm4hep::Vector3f(float((*p)->momentum().px()), float((*p)->momentum().py()), float((*p)->momentum().pz()) ));
        const HepMC::Polarization & polar = (*p)->polarization();
        mcp.setSpin               (edm4hep::Vector3f(polar.normal3d().x(), polar.normal3d().y(), polar.normal3d().z()) );
        int two[2] = {1, (*p)->flow(1)};
        mcp.setColorFlow          (edm4hep::Vector2i (two) );
    }
    // second loop for setting parents and daughters
    index = 0 ;
    for ( HepMC::GenEvent::particle_iterator p = evt->particles_begin(); p != evt->particles_end(); ++p ) {
        edm4hep::MCParticle pmc = event.m_mc_vec.at(index);
        index++;
        if ( (*p)->production_vertex() ) {
            for ( HepMC::GenVertex::particle_iterator mother = (*p)->production_vertex()-> particles_begin(HepMC::parents); mother != (*p)->production_vertex()-> particles_end(HepMC::parents); ++mother ) {
                pmc.addToParents( event.m_mc_vec.at( pmcid_lmcid.at((*mother)->barcode()) ) );
            }
        }
        if ( (*p)->end_vertex() ) {
            for ( HepMC::GenVertex::particle_iterator des =(*p)->end_vertex()-> particles_begin(HepMC::descendants); des != (*p)->end_vertex()-> particles_end(HepMC::descendants); ++des ) {
                pmc.addToDaughters( event.m_mc_vec.at( pmcid_lmcid.at((*des)->barcode()) ) );
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

bool HepMCRdr::configure_gentool(){
    ascii_in = new HepMC::IO_GenEvent(m_filename.value().c_str(),std::ios::in);

    m_processed_event=0;
    return true;
}

bool HepMCRdr::finish(){
    return true;
}

StatusCode
HepMCRdr::initialize() {
    StatusCode sc;
    if (not configure_gentool()) {
        error() << "failed to initialize." << endmsg;
        return StatusCode::FAILURE;
    }

    return sc;
}

StatusCode
HepMCRdr::finalize() {
    StatusCode sc;
    if (not finish()) {
        error() << "Failed to finalize." << endmsg;
        return StatusCode::FAILURE;
    }
    return sc;
}
