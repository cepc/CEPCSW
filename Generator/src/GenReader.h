#ifndef GenReader_h
#define GenReader_h 1

#include "GenEvent.h"
#include "IGenTool.h"

using namespace std;

class GenReader: virtual public IGenTool{

    public:
        ~GenReader();
        virtual bool configure_gentool()=0;               
        virtual bool mutate(MyHepMC::GenEvent& event)=0;    
        virtual bool finish()=0;
        virtual bool isEnd()=0;
};

#endif
