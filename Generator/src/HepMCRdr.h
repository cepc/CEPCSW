#ifndef HepMCRdr_h
#define HepMCRdr_h 1

#include "GenReader.h"
#include "GenEvent.h"

#include "HepMC/IO_GenEvent.h"//HepMC
#include "HepMC/GenEvent.h"


class HepMCRdr: public GenReader{

    public:
        HepMCRdr(string name);
        ~HepMCRdr();
        bool configure();               
        bool mutate(MyHepMC::GenEvent& event);    
        bool finish();
        bool isEnd();
    private:
        HepMC::IO_GenEvent *ascii_in;
        long m_total_event;
        long m_processed_event;
};

#endif

