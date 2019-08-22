#ifndef TEST_PLCIO_WRITE_ALG_H
#define TEST_PLCIO_WRITE_ALG_H

#include "FWCore/DataHandle.h"
#include "GaudiAlg/GaudiAlgorithm.h"

namespace plcio {
    class MCParticleCollection;
}

class PlcioWriteAlg : public GaudiAlgorithm
{
        friend class AlgFactory<PlcioWriteAlg>;

    public :

        PlcioWriteAlg(const std::string& name, ISvcLocator* svcLoc);

        virtual StatusCode initialize();
        virtual StatusCode execute();
        virtual StatusCode finalize();

    private :

        DataHandle<plcio::MCParticleCollection> m_hdl{"MCParticleCol", Gaudi::DataHandle::Writer, this};

};

#endif  // TEST_PLCIO_WRITE_ALG_H
