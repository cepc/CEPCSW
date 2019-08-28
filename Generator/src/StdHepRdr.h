#ifndef StdHepRdr_h
#define StdHepRdr_h 1

#include "GenReader.h"
#include "GenEvent.h"

#include "lcio.h"
#include "EVENT/LCIO.h"
#include "UTIL/LCStdHepRdrNew.h"


class StdHepRdr: public GenReader{

    public:
        StdHepRdr(string name);
        ~StdHepRdr();
        bool configure();               
        bool mutate(MyHepMC::GenEvent& event);    
        bool finish();
        bool isEnd();
    private:
        lcio::LCStdHepRdrNew* m_stdhep_rdr;
        long m_total_event;
        long m_processed_event;
};

#endif

