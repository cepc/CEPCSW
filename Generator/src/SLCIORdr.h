#ifndef SLCIORdr_h
#define SLCIORdr_h 1

#include "GenReader.h"
#include "GenEvent.h"

#include "lcio.h"
#include "LCIOSTLTypes.h"
#include "IOIMPL/LCFactory.h"
#include "EVENT/LCIO.h"
#include "EVENT/LCEvent.h"
#include "EVENT/LCCollection.h"
#include "IO/LCReader.h"


class SLCIORdr: public GenReader{

    public:
        SLCIORdr(string name);
        ~SLCIORdr();
        bool configure();               
        bool mutate(MyHepMC::GenEvent& event);    
        bool finish();
        bool isEnd();
    private:
        IO::LCReader* m_slcio_rdr;
        long m_total_event;
        long m_processed_event;
};

#endif

