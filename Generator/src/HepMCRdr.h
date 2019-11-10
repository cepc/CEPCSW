#ifndef HepMCRdr_h
#define HepMCRdr_h 1

#include "GaudiKernel/AlgTool.h"

#include "GenReader.h"
#include "GenEvent.h"

#include "HepMC/IO_GenEvent.h"//HepMC
#include "HepMC/GenEvent.h"


class HepMCRdr: public extends<AlgTool, GenReader> {

    public:
        using extends::extends;
        ~HepMCRdr();

        StatusCode initialize() override;
        StatusCode finalize() override;    

        bool configure_gentool();               
        bool mutate(MyHepMC::GenEvent& event);    
        bool finish();
        bool isEnd();
    private:
        HepMC::IO_GenEvent *ascii_in{nullptr};
        long m_total_event{-1};
        long m_processed_event{-1};

        // input file name
        Gaudi::Property<std::string> m_filename{this, "Input"};

};

#endif

