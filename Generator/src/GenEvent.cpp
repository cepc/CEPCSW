#include "GenEvent.h" 
#include "edm4hep/MCParticleCollection.h"//plico

// using namespace std;

namespace MyHepMC{

//GenEvent::GenEvent(){
GenEvent::GenEvent(edm4hep::MCParticleCollection& mcCol)
    : m_mc_vec(mcCol){

    m_event_id=-1;
    m_run_id=-1;
    m_time=-1;
    m_det_name="";

}
GenEvent::~GenEvent(){}

void GenEvent::SetEventHeader(long event_id_, long run_id_, float time_, const std::string& det_name_){
    m_event_id = event_id_;
    m_run_id = run_id_;
    m_time = time_;
    m_det_name = det_name_;
}
/*
void GenEvent::SetMCCollection(edm4hep::MCParticleCollection vec_){
m_mc_vec = vec_;
}
*/

edm4hep::MCParticleCollection& GenEvent::getMCVec(){
    return m_mc_vec;
}

long GenEvent::getID(){
return m_event_id;
}
long GenEvent::getRun() {return m_run_id;}
long GenEvent::getTime() {return m_time;}
std::string GenEvent::getName() {return m_det_name;}

void GenEvent::ReSet(){

    m_event_id=-1;
    m_run_id=-1;
    m_time=-1;
    m_det_name="";
    m_mc_vec.clear();
}
}
