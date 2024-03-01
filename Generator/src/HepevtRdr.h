#ifndef HepevtRdr_h
#define HepevtRdr_h 1

#include "GaudiKernel/AlgTool.h"

#include "GenReader.h"
#include "GenEvent.h"

#include "lcio.h"
#include "EVENT/LCIO.h"
#include "LCAscHepRdr.h"

class HepevtRdr: public extends<AlgTool, GenReader> {

    public:
        using extends::extends;
        ~HepevtRdr();

        StatusCode initialize() override;
        StatusCode finalize() override;    

        bool configure_gentool();               
        bool mutate(MyHepMC::GenEvent& event);    
        bool finish();
        bool isEnd();
    private:
        UTIL::LCAscHepRdr* m_hepevt_rdr = nullptr;
        long m_total_event = -1;
        long m_processed_event = -1;

        // input file name
        Gaudi::Property<std::string> m_filename{this, "Input"};
        Gaudi::Property<std::string> m_format{this, "Format"};
};

#endif

