#include "PlcioWriteAlg.h"
#include "plcio/MCParticleCollection.h"

DECLARE_COMPONENT(PlcioWriteAlg)

PlcioWriteAlg::PlcioWriteAlg(const std::string& name, ISvcLocator* svcLoc)
    : GaudiAlgorithm(name, svcLoc)
{
    declareProperty("MCParticleCol", m_hdl, "MCParticle collection (output)");
}

StatusCode PlcioWriteAlg::initialize()
{
    debug() << "begin initialize PlcioWriteAlg" << endmsg;
    return GaudiAlgorithm::initialize();
}

StatusCode PlcioWriteAlg::execute()
{
    debug() << "begin execute PlcioWriteAlg" << endmsg;

    auto mcCol = new plcio::MCParticleCollection;

    auto p1 = mcCol->create();
    auto p2 = mcCol->create();

    for ( int i = 0; i < 4; ++i ) {
        auto d = mcCol->create();
        d.addParent(p1);
        d.addParent(p2);
        p1.addDaughter(d);
        p2.addDaughter(d);
    }

    m_hdl.put(mcCol);

    return StatusCode::SUCCESS;
}

StatusCode PlcioWriteAlg::finalize()
{
    debug() << "begin finalize PlcioWriteAlg" << endmsg;
    return GaudiAlgorithm::finalize();
}
