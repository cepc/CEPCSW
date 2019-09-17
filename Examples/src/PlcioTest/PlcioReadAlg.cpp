#include "PlcioReadAlg.h"
#include "plcio/MCParticleCollection.h"

DECLARE_COMPONENT(PlcioReadAlg)

PlcioReadAlg::PlcioReadAlg(const std::string& name, ISvcLocator* svcLoc)
    : GaudiAlgorithm(name, svcLoc)
{
    declareProperty("MCParticleCol", m_hdl, "MCParticle collection (input)");
}

StatusCode PlcioReadAlg::initialize()
{
    debug() << "begin initialize PlcioReadAlg" << endmsg;
    return GaudiAlgorithm::initialize();
}

StatusCode PlcioReadAlg::execute()
{
    debug() << "begin execute PlcioReadAlg" << endmsg;
    auto mcCol = m_hdl.get();

//    debug() << "testing loop..." <<endmsg;
    for ( auto p : *mcCol ) {
        debug() << p.getObjectID().index << " : [";
        for ( auto it = p.daughters_begin(), end = p.daughters_end(); it != end; ++it ) {
            debug() << " " << it->getObjectID().index;
        }
        debug() << " ]; ";
    }
    debug() << endmsg;

//    debug() << "end execute PlcioReadAlg" << endmsg;
    return StatusCode::SUCCESS;
}

StatusCode PlcioReadAlg::finalize()
{
    debug() << "begin finalize PlcioReadAlg" << endmsg;
    return GaudiAlgorithm::finalize();
}
