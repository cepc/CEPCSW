#ifndef GenWriter_h
#define GenWriter_h 1

#include "GenEvent.h"
#include "IGenTool.h"

#include "podio/EventStore.h" //podio
#include "podio/ROOTWriter.h"


#include "plcio/MCParticleCollection.h"//plico
#include "plcio/EventHeaderCollection.h"



using namespace std;

class GenWriter: public IGenTool{

    public:
        GenWriter(string name);
        ~GenWriter();
        bool configure() override;               
        bool mutate(MyHepMC::GenEvent& event) override;    
        bool finish() override;
    private:
        string m_output_name;
        podio::EventStore* store ;
        podio::ROOTWriter* writer;
        plcio::EventHeaderCollection*  ehc  ;
        plcio::MCParticleCollection*   mcc  ;

};

#endif
