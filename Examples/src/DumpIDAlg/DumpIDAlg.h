#ifndef DumpIDAlg_h
#define DumpIDAlg_h

#include "FWCore/DataHandle.h"
#include "GaudiAlg/GaudiAlgorithm.h"

#include "DetInterface/IGeoSvc.h"

#include "DD4hep/Detector.h"

namespace edm4hep {
    class EventHeaderCollection;
    class MCParticleCollection;
    class SimCalorimeterHitCollection;
    class CaloHitContributionCollection;
}

class DumpIDAlg: public GaudiAlgorithm
{
public:

    DumpIDAlg(const std::string& name, ISvcLocator* svcLoc);

    virtual StatusCode initialize();
    virtual StatusCode execute();
    virtual StatusCode finalize();

private:
    SmartIF<IGeoSvc> m_geosvc;
    dd4hep::Detector* m_dd4hep_geo;
    dd4hep::DDSegmentation::BitFieldCoder* m_decoder;

    DataHandle<edm4hep::SimCalorimeterHitCollection> m_EcalBarrelCol{"EcalBarrelCollection", 
            Gaudi::DataHandle::Reader, this};

};


#endif
