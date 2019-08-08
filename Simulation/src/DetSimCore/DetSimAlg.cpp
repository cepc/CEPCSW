#include "DetSimAlg.h"

DECLARE_COMPONENT(DetSimAlg)

DetSimAlg::DetSimAlg(const std::string& name, ISvcLocator* pSvcLocator)
: Algorithm(name, pSvcLocator) {

}

StatusCode
DetSimAlg::initialize() {
    StatusCode sc;

    info() << "Initialize DetSimAlg... " << endmsg;

    return sc;
}

StatusCode
DetSimAlg::execute() {
    StatusCode sc;
    return sc;
}

StatusCode
DetSimAlg::finalize() {
    StatusCode sc;
    return sc;
}



