#include "GenPrinter.h"
#include "GenEvent.h"

DECLARE_COMPONENT(GenPrinter)

bool GenPrinter::mutate(MyHepMC::GenEvent& event){
    std::cout << "print mc info for event "<< event.getID() << ", mc size ="<< event.m_mc_vec.size() <<  std::endl;
    for ( int i =0; i < event.m_mc_vec.size(); i++ ) {
    auto p = event.m_mc_vec.at(i); 
    std::cout<< "PDG      :"<< p.getPDG               ()<<std::endl 
    << "id                :"<< p.id                   ()<<std::endl 
    << "ID                :"<< p.getObjectID().index    <<std::endl 
    << "GeneratorStatus   :"<< p.getGeneratorStatus   ()<<std::endl 
    << "SimulatorStatus   :"<< p.getSimulatorStatus   ()<<std::endl 
    << "Charge            :"<< p.getCharge            ()<<std::endl 
    << "Time              :"<< p.getTime              ()<<std::endl 
    << "Mass              :"<< p.getMass              ()<<std::endl 
    << "Vertex            :"<< p.getVertex            ()<<std::endl 
    << "Endpoint          :"<< p.getEndpoint          ()<<std::endl 
    << "Momentum          :"<< p.getMomentum          ()<<std::endl 
    << "MomentumAtEndpoint:"<< p.getMomentumAtEndpoint()<<std::endl 
    << "Spin              :"<< p.getSpin              ()<<std::endl 
    << "ColorFlow         :"<< p.getColorFlow         ()<<std::endl 
    << "Parent size       :"<< p.parents_size         ()<<std::endl 
    << "Daughter size     :"<< p.daughters_size       ()<<std::endl; 
    //for(unsigned int j=0; j<p.parents_size(); j++) std::cout << " for parent: "<< j << ",PDG="<< p.getParents(j).getPDG() << ",id=:"<< p.getParents(j).id()<<std::endl;
    for (auto it = p.parents_begin(), end = p.parents_end(); it != end ; ++it ) std::cout << " for parent ,PDG="<< it->getPDG() << ",id=:"<< it->getObjectID().index<<std::endl;
    }
    return true;
}

bool GenPrinter::configure_gentool(){
    return true;
}

bool GenPrinter::finish(){
    return true;
}

StatusCode
GenPrinter::initialize() {
    StatusCode sc;
    if (not configure_gentool()) {
        error() << "failed to initialize." << endmsg;
        return StatusCode::FAILURE;
    }

    return sc;
}

StatusCode
GenPrinter::finalize() {
    StatusCode sc;
    if (not finish()) {
        error() << "Failed to finalize." << endmsg;
        return StatusCode::FAILURE;
    }
    return sc;
}

