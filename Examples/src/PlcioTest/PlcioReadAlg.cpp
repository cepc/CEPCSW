#include "PlcioReadAlg.h"
#include "plcio/EventHeaderCollection.h"
#include "plcio/MCParticleCollection.h"

DECLARE_COMPONENT(PlcioReadAlg)

PlcioReadAlg::PlcioReadAlg(const std::string& name, ISvcLocator* svcLoc)
    : GaudiAlgorithm(name, svcLoc)
{
    declareProperty("HeaderCol", m_headerCol);
    declareProperty("InputCol", m_mcParCol, "MCParticle collection (input)");
}

StatusCode PlcioReadAlg::initialize()
{
    debug() << "begin initialize PlcioReadAlg" << endmsg;
    return GaudiAlgorithm::initialize();
}

StatusCode PlcioReadAlg::execute()
{
    debug() << "begin execute PlcioReadAlg" << endmsg;

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

StatusCode PlcioReadAlg::finalize()
{
    debug() << "begin finalize PlcioReadAlg" << endmsg;
    return GaudiAlgorithm::finalize();
}
