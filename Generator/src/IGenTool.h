#ifndef IGenTool_h
#define IGenTool_h 1

/*
 * IGenTool is used to mutate an event. Several tools could be used 
 * together to create a full event.
 *
 * ChangeLog:
 * - 2019.11.10, Tao Lin, make the IGenTool a Gaudi tool.
 */

#include "GaudiKernel/IAlgTool.h"
#include "GenEvent.h"


class IGenTool: virtual public IAlgTool  {
    public:
        virtual bool mutate(MyHepMC::GenEvent& event)=0;
        virtual bool finish()=0;
        virtual bool configure_gentool()=0;
        virtual ~IGenTool();
};

#endif
