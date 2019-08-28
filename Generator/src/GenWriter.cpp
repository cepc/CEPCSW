#include "GenWriter.h"
#include "GenEvent.h"

#include "podio/EventStore.h" //podio
#include "podio/ROOTWriter.h"


#include "plcio/MCParticleCollection.h"//plico
#include "plcio/EventHeaderCollection.h"


GenWriter::GenWriter(string name){
    m_output_name = name;
    store  = new podio::EventStore();
    writer = new podio::ROOTWriter(m_output_name, store);
    ehc   =  &store->create<plcio::EventHeaderCollection>("EvtHeaders");
    mcc   =  &store->create<plcio::MCParticleCollection>("MCParticles");

    writer->registerForWrite("EvtHeaders");
    writer->registerForWrite("MCParticles");
}

GenWriter::~GenWriter(){
}

bool GenWriter::mutate(MyHepMC::GenEvent& event){
    std::cout << "write mc info for event "<< event.getID() << ", mc size ="<< event.m_mc_vec.size() <<  std::endl;
    mcc=&event.m_mc_vec;
    auto header = plcio::EventHeader(event.getID(), event.getRun(), event.getTime(), event.getName());
    ehc->push_back(header);
    writer->writeEvent();
    store->clearCollections();
    return true;
}

bool GenWriter::configure(){
return true;
}

bool GenWriter::finish(){
    writer->finish();
    std::cout<<"Saved root "<<m_output_name<<std::endl;
    return true;
}
