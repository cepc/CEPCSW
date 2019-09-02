#include "PlcioWriteAlg.h"
#include "plcio/EventHeaderCollection.h"
#include "plcio/MCParticleCollection.h"

DECLARE_COMPONENT(PlcioWriteAlg)

PlcioWriteAlg::PlcioWriteAlg(const std::string& name, ISvcLocator* svcLoc)
    : GaudiAlgorithm(name, svcLoc)
{
    declareProperty("HeaderCol", m_headerCol);
    declareProperty("OutputCol", m_mcParCol, "MCParticle collection (output)");
}

StatusCode PlcioWriteAlg::initialize()
{
    debug() << "begin initialize PlcioWriteAlg" << endmsg;
    return GaudiAlgorithm::initialize();
}

StatusCode PlcioWriteAlg::execute()
{
    debug() << "begin execute PlcioWriteAlg" << endmsg;

    static int evtNo = 0;

    auto headers = m_headerCol.createAndPut();
    auto header = headers->create();
    header->setRunNumber(-999);
    header->setEventNumber(evtNo++);
    header->setDetectorName("TEST");

    //auto mcCol = new plcio::MCParticleCollection;
    //m_mcParCol.put(mcCol);
    auto mcCol = m_mcParCol.createAndPut();

    auto p1 = mcCol->create();
    auto p2 = mcCol->create();

    for ( int i = 0; i < 4; ++i ) {
        auto d = mcCol->create();
        d.addParent(p1);
        d.addParent(p2);
        p1.addDaughter(d);
        p2.addDaughter(d);
    }

    return StatusCode::SUCCESS;
}

StatusCode PlcioWriteAlg::finalize()
{
    debug() << "begin finalize PlcioWriteAlg" << endmsg;
    return GaudiAlgorithm::finalize();
}
