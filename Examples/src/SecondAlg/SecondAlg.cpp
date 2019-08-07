#include "SecondAlg.h"

DECLARE_COMPONENT(SecondAlg)

SecondAlg::SecondAlg(const std::string& name, ISvcLocator* pSvcLocator)
: Algorithm(name, pSvcLocator) {

}

StatusCode
SecondAlg::initialize() {
    StatusCode sc;

    info() << "Retrieving the FirstSvc... " << endmsg;

    m_firstsvc = service("FirstSvc");

    return sc;
}

StatusCode
SecondAlg::execute() {
    StatusCode sc;

    m_firstsvc->shoot();

    return sc;
}

StatusCode
SecondAlg::finalize() {
    StatusCode sc;
    return sc;
}


