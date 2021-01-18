#ifndef TEST_EDM4HEP_READ_DC_ALG_H
#define TEST_EDM4HEP_READ_DC_ALG_H

#include "k4FWCore/DataHandle.h"
#include "GaudiAlg/GaudiAlgorithm.h"

namespace edm4hep {
    class EventHeaderCollection;
    class MCParticleCollection;
    class SimTrackerHitCollection;
}

class Edm4hepReadDCAlg : public GaudiAlgorithm
{

    public :

        Edm4hepReadDCAlg(const std::string& name, ISvcLocator* svcLoc);

        virtual StatusCode initialize();
        virtual StatusCode execute();
        virtual StatusCode finalize();

    private :

        // DataHandle<edm4hep::EventHeaderCollection> m_headerCol{"EventHeader", Gaudi::DataHandle::Reader, this};
        DataHandle<edm4hep::MCParticleCollection> m_mcParCol{"MCParticle", Gaudi::DataHandle::Reader, this};
        DataHandle<edm4hep::SimTrackerHitCollection> m_dcCol{"DriftChamberHitsCollection", 
                Gaudi::DataHandle::Reader, this};

};

#endif
