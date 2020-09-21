#include "DummyDedxSimTool.h"

#include "G4Step.hh"

DECLARE_COMPONENT(DummyDedxSimTool);

StatusCode DummyDedxSimTool::initialize() {
    StatusCode sc;

    return sc;
}

StatusCode DummyDedxSimTool::finalize() {
    StatusCode sc;

    return sc;
}

double DummyDedxSimTool::dedx(const G4Step* aStep) {
    double result = aStep->GetTotalEnergyDeposit();


    return result;
}
