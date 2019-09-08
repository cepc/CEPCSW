#ifndef G4PrimaryCnvTool_h
#define G4PrimaryCnvTool_h

#include "GaudiKernel/AlgTool.h"
#include "DetSimInterface/IG4PrimaryCnvTool.h"
#include "FWCore/DataHandle.h"

#include "plcio/EventHeaderCollection.h"
#include "plcio/MCParticleCollection.h"

class G4PrimaryCnvTool: public extends<AlgTool, IG4PrimaryCnvTool> {
public:

    using extends::extends;

    bool mutate(G4Event* anEvent) override;

private:
    DataHandle<plcio::MCParticleCollection> m_mcParCol{"MCParticleCol", Gaudi::DataHandle::Reader, this};

};

#endif
