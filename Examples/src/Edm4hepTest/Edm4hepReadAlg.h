#ifndef TEST_EDM4HEP_WRITE_ALG_H
#define TEST_EDM4HEP_WRITE_ALG_H

#include "FWCore/DataHandle.h"
#include "GaudiAlg/GaudiAlgorithm.h"

namespace edm4hep {
    class EventHeaderCollection;
    class MCParticleCollection;
}

class Edm4hepReadAlg : public GaudiAlgorithm
{

    public :

        Edm4hepReadAlg(const std::string& name, ISvcLocator* svcLoc);

        virtual StatusCode initialize();
        virtual StatusCode execute();
        virtual StatusCode finalize();

    private :

        DataHandle<edm4hep::EventHeaderCollection> m_headerCol{"EventHeader", Gaudi::DataHandle::Reader, this};
        DataHandle<edm4hep::MCParticleCollection> m_mcParCol{"MCParticle", Gaudi::DataHandle::Reader, this};

};

#endif  // TEST_EDM4HEP_WRITE_ALG_H
