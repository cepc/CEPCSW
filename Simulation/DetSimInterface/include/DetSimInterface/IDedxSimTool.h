#ifndef IDedxSimTool_h
#define IDedxSimTool_h

/*
 * Description:
 *   IDedxSimTool is used to give a dE/dx value during simulation.
 *
 * The interface:
 *   * dedx: predict the dE/dx according to Geant4 Step
 *
 * Author: Tao Lin <lintao@ihep.ac.cn>
 */


#include "GaudiKernel/IAlgTool.h"

class G4Step;
namespace edm4hep{
    class MCParticle;
}

class IDedxSimTool: virtual public IAlgTool {
public:

    DeclareInterfaceID(IDedxSimTool, 0, 1);
    virtual ~IDedxSimTool() {}

    virtual double dedx(const G4Step* aStep) = 0;
    virtual double dedx(const edm4hep::MCParticle& mc) = 0;
    virtual double dndx(double betagamma) = 0;

};

#endif
