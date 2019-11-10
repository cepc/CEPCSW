#ifndef SLCIORdr_h
#define SLCIORdr_h 1

#include "GaudiKernel/AlgTool.h"

#include "GenReader.h"
#include "GenEvent.h"

#include "lcio.h"
#include "LCIOSTLTypes.h"
#include "IOIMPL/LCFactory.h"
#include "EVENT/LCIO.h"
#include "EVENT/LCEvent.h"
#include "EVENT/LCCollection.h"
#include "IO/LCReader.h"


class SLCIORdr: public extends<AlgTool, GenReader> {

    public:
        using extends::extends;

        ~SLCIORdr();

        // Overriding initialize and finalize
        StatusCode initialize() override;
        StatusCode finalize() override;    

        bool configure_gentool() override;
        bool mutate(MyHepMC::GenEvent& event) override;    
        bool finish() override;
        bool isEnd() override;
    private:
        IO::LCReader* m_slcio_rdr{nullptr};
        long m_total_event{-1};
        long m_processed_event{-1};

        // input file name
        Gaudi::Property<std::string> m_filename{this, "Input"};

};

#endif

