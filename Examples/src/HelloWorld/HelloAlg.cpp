#include "HelloAlg.h"

DECLARE_COMPONENT(HelloAlg)

HelloAlg::HelloAlg(const std::string& name, ISvcLocator* pSvcLocator)
: Algorithm(name, pSvcLocator) {

}

StatusCode
HelloAlg::initialize() {
    StatusCode sc;

    info() << "MyInt: " << m_int.value() << endmsg;

    return sc;
}

StatusCode
HelloAlg::execute() {
    StatusCode sc;
    return sc;
}

StatusCode
HelloAlg::finalize() {
    StatusCode sc;
    return sc;
}


