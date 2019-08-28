#ifndef HepevtRdr_h
#define HepevtRdr_h 1

#include "GenReader.h"
#include "GenEvent.h"

#include "lcio.h"
#include "EVENT/LCIO.h"
#include "LCAscHepRdr.h"

class HepevtRdr: public GenReader{

    public:
        HepevtRdr(string name);
        ~HepevtRdr();
        bool configure();               
        bool mutate(MyHepMC::GenEvent& event);    
        bool finish();
        bool isEnd();
    private:
        UTIL::LCAscHepRdr* m_hepevt_rdr;
        long m_total_event;
        long m_processed_event;
};

#endif

