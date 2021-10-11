#include "GtBeamBackgroundTool.h"

DECLARE_COMPONENT(GtBeamBackgroundTool)

StatusCode GtBeamBackgroundTool::initialize() {
    // check the consistency of the properties

    // create the instances of the background parsers

    return StatusCode::SUCCESS;
}

StatusCode GtBeamBackgroundTool::finalize() {
    return StatusCode::SUCCESS;
}


bool GtBeamBackgroundTool::mutate(MyHepMC::GenEvent& event) {
    return true;
}

bool GtBeamBackgroundTool::finish() {
    return true;
}

bool GtBeamBackgroundTool::configure_gentool() {

    return true;
}
