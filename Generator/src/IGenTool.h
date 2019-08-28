#ifndef IGenTool_h
#define IGenTool_h 1

#include "GenEvent.h"


class IGenTool {
    public:
        virtual bool mutate(MyHepMC::GenEvent& event)=0;
        virtual bool finish()=0;
        virtual bool configure()=0;
        virtual ~IGenTool();
};

#endif
