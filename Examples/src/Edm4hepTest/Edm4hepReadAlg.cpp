#include "Edm4hepReadAlg.h"
#include "edm4hep/EventHeaderCollection.h"
#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/SimCalorimeterHitCollection.h"
#include "edm4hep/CaloHitContributionCollection.h"

DECLARE_COMPONENT(Edm4hepReadAlg)

Edm4hepReadAlg::Edm4hepReadAlg(const std::string& name, ISvcLocator* svcLoc)
    : GaudiAlgorithm(name, svcLoc)
{
    declareProperty("HeaderCol", m_headerCol);
    declareProperty("MCParticleCol", m_mcParCol, "MCParticle collection (input)");
    declareProperty("SimCalorimeterHitCol", m_calorimeterCol, "MCParticle collection (input)");
}

StatusCode Edm4hepReadAlg::initialize()
{
    debug() << "begin initialize Edm4hepReadAlg" << endmsg;
    return GaudiAlgorithm::initialize();
}

StatusCode Edm4hepReadAlg::execute()
{
    debug() << "begin execute Edm4hepReadAlg" << endmsg;

    auto mcCol = m_mcParCol.get();
    for ( auto p : *mcCol ) {
        info() << p.getObjectID().index << " : [";
        for ( auto it = p.daughters_begin(), end = p.daughters_end(); it != end; ++it ) {
            info() << " " << it->getObjectID().index;
        }
        info() << " ]; ";
    }
    info() << "}" << endmsg;

    auto caloCol = m_calorimeterCol.get();
    for (auto calohit : *caloCol) {
        unsigned int contrib_size = calohit.contributions_size();
        info() << " contributions_size: " 
               << contrib_size
               << endmsg;
        for (unsigned int i = 0; i < contrib_size; ++i) {
            auto contrib = calohit.getContributions(i);
            auto primary_particle = contrib.getParticle();

            info() << " - #" << i << ": "
                   << " track with "
                   << " PDG: " << contrib.getPDG() // current track
                   << ". "
                   << " primary track with "
                   << " PDG: " << primary_particle.getPDG()
                   << endmsg;


        }

    }

    return StatusCode::SUCCESS;
}

StatusCode Edm4hepReadAlg::finalize()
{
    debug() << "begin finalize Edm4hepReadAlg" << endmsg;
    return GaudiAlgorithm::finalize();
}
