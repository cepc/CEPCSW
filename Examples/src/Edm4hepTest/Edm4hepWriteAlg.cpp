#include "Edm4hepWriteAlg.h"
#include "edm4hep/EventHeaderCollection.h"
#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/SimCalorimeterHitCollection.h"
#include "edm4hep/CaloHitContributionCollection.h"

DECLARE_COMPONENT(Edm4hepWriteAlg)

Edm4hepWriteAlg::Edm4hepWriteAlg(const std::string& name, ISvcLocator* svcLoc)
    : GaudiAlgorithm(name, svcLoc)
{
    declareProperty("HeaderOut", m_headerCol);
    declareProperty("MCParticleOut", m_mcParCol, "MCParticle collection (output)");
    declareProperty("SimCalorimeterHitOut", m_simCaloHitCol, "SimCalorimeterHit collection (output)");
    declareProperty("CaloHitContributionOut", m_caloHitContCol, "CaloHitContribution collection (output)");
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
    auto simCaloCol = m_simCaloHitCol.createAndPut();
    auto caloHitContCol = m_caloHitContCol.createAndPut();

    auto p1 = mcCol->create();
    auto p2 = mcCol->create();

    for ( int i = 0; i < 4; ++i ) {
        auto d = mcCol->create();
        d.addToParents(p1);
        d.addToParents(p2);
        p1.addToDaughters(d);
        p2.addToDaughters(d);

        auto hit = simCaloCol->create();
        for ( int j = 0; j < i; ++j ) {
            auto cont = caloHitContCol->create();
            cont.setParticle(mcCol->at(j));
            hit.addToContributions(cont);
        }
    }

    return StatusCode::SUCCESS;
}

StatusCode Edm4hepWriteAlg::finalize()
{
    debug() << "begin finalize Edm4hepWriteAlg" << endmsg;
    return GaudiAlgorithm::finalize();
}
