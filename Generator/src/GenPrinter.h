#ifndef GenPrinter_h
#define GenPrinter_h 1

#include "GenEvent.h"
#include "IGenTool.h"

using namespace std;

class GenPrinter: public IGenTool{

    public:
        GenPrinter(string name);
        ~GenPrinter();
        bool configure() override;               
        bool mutate(MyHepMC::GenEvent& event) override;    
        bool finish() override;
};

#endif
