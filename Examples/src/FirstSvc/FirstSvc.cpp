#include "FirstSvc.h"

DECLARE_COMPONENT(FirstSvc)

StatusCode
FirstSvc::initialize() {
    StatusCode sc = Service::initialize();

    return sc;
}

StatusCode
FirstSvc::finalize() {
    // clear or reset

    // 
    StatusCode sc = Service::finalize();
    return sc;
}

void
FirstSvc::shoot() {
    msg() << "shoot it." << endmsg;

}
