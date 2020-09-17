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
    double result = 0.0;


    return result;
}
