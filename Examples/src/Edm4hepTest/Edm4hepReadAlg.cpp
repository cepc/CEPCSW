#include "Edm4hepReadAlg.h"
#include "edm4hep/EventHeaderCollection.h"
#include "edm4hep/MCParticleCollection.h"

DECLARE_COMPONENT(Edm4hepReadAlg)

Edm4hepReadAlg::Edm4hepReadAlg(const std::string& name, ISvcLocator* svcLoc)
    : GaudiAlgorithm(name, svcLoc)
{
    declareProperty("HeaderCol", m_headerCol);
    declareProperty("InputCol", m_mcParCol, "MCParticle collection (input)");
}

StatusCode Edm4hepReadAlg::initialize()
{
    debug() << "begin initialize Edm4hepReadAlg" << endmsg;
    return GaudiAlgorithm::initialize();
}

StatusCode Edm4hepReadAlg::execute()
{
    debug() << "begin execute Edm4hepReadAlg" << endmsg;

    auto headers = m_headerCol.get();
    auto header = headers->at(0);
    auto mcCol = m_mcParCol.get();

    info() << "Run " << header.getRunNumber() << " Event " << header.getEventNumber() << " { ";
    for ( auto p : *mcCol ) {
        info() << p.getObjectID().index << " : [";
        for ( auto it = p.daughters_begin(), end = p.daughters_end(); it != end; ++it ) {
            info() << " " << it->getObjectID().index;
        }
        info() << " ]; ";
    }
    info() << "}" << endmsg;

    return StatusCode::SUCCESS;
}

StatusCode Edm4hepReadAlg::finalize()
{
    debug() << "begin finalize Edm4hepReadAlg" << endmsg;
    return GaudiAlgorithm::finalize();
}
