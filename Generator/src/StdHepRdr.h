#ifndef StdHepRdr_h
#define StdHepRdr_h 1

#include "GaudiKernel/AlgTool.h"

#include "GenReader.h"
#include "GenEvent.h"

#include "lcio.h"
#include "EVENT/LCIO.h"
#include "UTIL/LCStdHepRdrNew.h"


class StdHepRdr: public extends<AlgTool, GenReader> {

public:

    using extends::extends;

    ~StdHepRdr();

    // Overriding initialize and finalize
    StatusCode initialize() override;
    StatusCode finalize() override;    

    bool configure_gentool() override;               
    bool mutate(MyHepMC::GenEvent& event) override;    
    bool finish() override;
    bool isEnd() override;
private:
    lcio::LCStdHepRdrNew* m_stdhep_rdr{nullptr};
    long m_total_event{-1};
    long m_processed_event{-1};

    // input file name
    Gaudi::Property<std::string> m_filename{this, "Input"};

};

#endif

