#ifndef G4PrimaryCnvTool_h
#define G4PrimaryCnvTool_h

#include "GaudiKernel/AlgTool.h"
#include "DetSimInterface/IG4PrimaryCnvTool.h"
#include "k4FWCore/DataHandle.h"

#include "edm4hep/EventHeaderCollection.h"
#include "edm4hep/MCParticleCollection.h"

class G4PrimaryCnvTool: public extends<AlgTool, IG4PrimaryCnvTool> {
public:

    using extends::extends;

    bool mutate(G4Event* anEvent) override;

private:
    DataHandle<edm4hep::MCParticleCollection> m_mcParCol{"MCParticleGen", Gaudi::DataHandle::Reader, this};

    Gaudi::Property<double> m_chargedgeantino_mass{this, "ChargedGeantinoMass"};
};

#endif
