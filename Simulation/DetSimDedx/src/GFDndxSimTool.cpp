#include "GFDndxSimTool.h"

#include "G4Step.hh"

DECLARE_COMPONENT(GFDndxSimTool);

StatusCode GFDndxSimTool::initialize() {
    StatusCode sc;

    return sc;
}

StatusCode GFDndxSimTool::finalize() {
    StatusCode sc;

    return sc;
}

double GFDndxSimTool::dedx(const G4Step* aStep) {
    double result = aStep->GetTotalEnergyDeposit();
    return result;
}
double GFDndxSimTool::dedx(const edm4hep::MCParticle& mc) {
    return -1;
}
double GFDndxSimTool::dndx(double betagamma) {
    return -999;
}
