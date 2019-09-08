#include "G4PrimaryCnvTool.h"

#include "G4Event.hh"

DECLARE_COMPONENT(G4PrimaryCnvTool)

bool G4PrimaryCnvTool::mutate(G4Event* anEvent) {

    auto mcCol = m_mcParCol.get();
    info() << "Start a new event: " << endmsg;
    for ( auto p : *mcCol ) {
        info() << p.getObjectID().index << " : [";
        for ( auto it = p.daughters_begin(), end = p.daughters_end(); it != end; ++it ) {
            info() << " " << it->getObjectID().index;
        }
        info() << " ]; " << endmsg;
    }

    return true;
}
