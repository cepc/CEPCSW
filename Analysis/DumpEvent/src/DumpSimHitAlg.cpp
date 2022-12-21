/*
 * Description:
 *   Dump the simulated information.
 *
 * Author:
 *   Tao Lin <lintao AT ihep.ac.cn>
 */

#include "k4FWCore/DataHandle.h"
#include "GaudiKernel/Algorithm.h"

#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/SimTrackerHitCollection.h"
#include "edm4hep/SimCalorimeterHitCollection.h"
#include "edm4hep/CaloHitContributionCollection.h"

#include "GaudiKernel/NTuple.h"

class DumpSimHitAlg: public Algorithm {
public:

    DumpSimHitAlg(const std::string& name, ISvcLocator* pSvcLocator);

    // Three mandatory member functions of any algorithm
    StatusCode initialize() override;
    StatusCode execute() override;
    StatusCode finalize() override;

private:
    // - collection MCParticleG4: the simulated particles in Geant4
    DataHandle<edm4hep::MCParticleCollection> m_mcParCol{"MCParticle", 
            Gaudi::DataHandle::Reader, this};

    // Dedicated collections for CEPC
    DataHandle<edm4hep::SimTrackerHitCollection> m_VXDCol{"VXDCollection", 
            Gaudi::DataHandle::Reader, this};

};

DECLARE_COMPONENT( DumpSimHitAlg )

DumpSimHitAlg::DumpSimHitAlg(const std::string& name, ISvcLocator* pSvcLocator) 
: Algorithm(name, pSvcLocator) {

}

StatusCode DumpSimHitAlg::initialize() {
    return StatusCode::SUCCESS;
}

StatusCode DumpSimHitAlg::execute() {
    auto mcCol = m_mcParCol.get();

    for (auto particle: *mcCol) {
        info() << "mc particle -> "
               << " (ID: " << particle.getObjectID().index << ") "
               << " (simulator status: " << particle.getSimulatorStatus() << ") "
               << endmsg;
    }

    auto vxdCol = m_VXDCol.get();

    for (auto hit: *vxdCol) {
        auto mcparticle = hit.getMCParticle();

        if (mcparticle.getGeneratorStatus() != 1) {
            error() << "Found generator status is not 1 for hit. " << endmsg;
        }

        info() << "hit -> "
               << " mcparticle ("
               << " ID: " << mcparticle.getObjectID().index << "; "
               << " generator status: " << mcparticle.getGeneratorStatus() << "; "
               << " simulator status: " << mcparticle.getSimulatorStatus() << ") "
               << endmsg;
    }

    return StatusCode::SUCCESS;
}

StatusCode DumpSimHitAlg::finalize() {
    return StatusCode::SUCCESS;
}
