#ifndef GenPrinter_h
#define GenPrinter_h 1

#include <GaudiKernel/AlgTool.h>
#include "GenEvent.h"
#include "IGenTool.h"

using namespace std;

class GenPrinter: public extends<AlgTool, IGenTool> {
public:
    using extends::extends;

    // Overriding initialize and finalize
    StatusCode initialize() override;
    StatusCode finalize() override;

public:
    bool configure_gentool() override;               
    bool mutate(MyHepMC::GenEvent& event) override;    
    bool finish() override;
};

#endif
