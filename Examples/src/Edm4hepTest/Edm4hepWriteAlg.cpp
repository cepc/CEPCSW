#include "Edm4hepWriteAlg.h"
#include "edm4hep/EventHeaderCollection.h"
#include "edm4hep/MCParticleCollection.h"

DECLARE_COMPONENT(Edm4hepWriteAlg)

Edm4hepWriteAlg::Edm4hepWriteAlg(const std::string& name, ISvcLocator* svcLoc)
    : GaudiAlgorithm(name, svcLoc)
{
    declareProperty("HeaderCol", m_headerCol);
    declareProperty("OutputCol", m_mcParCol, "MCParticle collection (output)");
}

StatusCode Edm4hepWriteAlg::initialize()
{
    debug() << "begin initialize Edm4hepWriteAlg" << endmsg;
    return GaudiAlgorithm::initialize();
}

StatusCode Edm4hepWriteAlg::execute()
{
    debug() << "begin execute Edm4hepWriteAlg" << endmsg;

    static int evtNo = 0;

    auto headers = m_headerCol.createAndPut();
    auto header = headers->create();
    header.setRunNumber(-999);
    header.setEventNumber(evtNo++);
    // header.setDetectorName("TEST");

    //auto mcCol = new edm4hep::MCParticleCollection;
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

StatusCode Edm4hepWriteAlg::finalize()
{
    debug() << "begin finalize Edm4hepWriteAlg" << endmsg;
    return GaudiAlgorithm::finalize();
}
