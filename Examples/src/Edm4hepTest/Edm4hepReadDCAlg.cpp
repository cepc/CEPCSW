#include "Edm4hepReadDCAlg.h"
#include "edm4hep/EventHeaderCollection.h"
#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/SimTrackerHitCollection.h"

DECLARE_COMPONENT(Edm4hepReadDCAlg)

Edm4hepReadDCAlg::Edm4hepReadDCAlg(const std::string& name, ISvcLocator* svcLoc)
    : GaudiAlgorithm(name, svcLoc)
{
    declareProperty("MCParticleCol", m_mcParCol, "MCParticle collection (input)");
    declareProperty("DCHitCol", m_dcCol, "Drift Chamber collections (input)");
}

StatusCode Edm4hepReadDCAlg::initialize()
{
    debug() << "begin initialize Edm4hepReadDCAlg" << endmsg;
    return GaudiAlgorithm::initialize();
}

StatusCode Edm4hepReadDCAlg::execute()
{
    debug() << "begin execute Edm4hepReadDCAlg" << endmsg;

    auto mcCol = m_mcParCol.get();
    for ( auto p : *mcCol ) {
        info() << p.getObjectID().index << " : [";
        for ( auto it = p.daughters_begin(), end = p.daughters_end(); it != end; ++it ) {
            info() << " " << it->getObjectID().index;
        }
        info() << " ]; ";
    }
    info() << "}" << endmsg;

    auto trkCol = m_dcCol.get();
    for (auto trkhit : *trkCol) {
        auto position = trkhit.getPosition();
        auto momentum = trkhit.getMomentum();
        auto primary_particle = trkhit.getMCParticle();

        info() << " cellID: " << trkhit.getCellID()
               << " edep: " << trkhit.getEDep()
               << " time: " << trkhit.getTime()
               << " length: " << trkhit.getPathLength()
               << " quality: " << trkhit.getQuality()
               << " position: ("
               << position[0] << ", " << position[1] << ", " << position[2]
               << ")"
               << " momentum: ("
               << momentum[0] << ", " << momentum[1] << ", " << momentum[2]
               << ")"
               << " primary track: "
               << primary_particle.getPDG()
               << endmsg;
            

        // unsigned int contrib_size = calohit.contributions_size();
        // info() << " contributions_size: " 
        //        << contrib_size
        //        << endmsg;
        // for (unsigned int i = 0; i < contrib_size; ++i) {
        //     auto contrib = calohit.getContributions(i);
        //     auto primary_particle = contrib.getParticle();

        //     info() << " - #" << i << ": "
        //            << " track with "
        //            << " PDG: " << contrib.getPDG() // current track
        //            << ". "
        //            << " primary track with "
        //            << " PDG: " << primary_particle.getPDG()
        //            << endmsg;


        // }

    }

    return StatusCode::SUCCESS;
}

StatusCode Edm4hepReadDCAlg::finalize()
{
    debug() << "begin finalize Edm4hepReadDCAlg" << endmsg;
    return GaudiAlgorithm::finalize();
}
