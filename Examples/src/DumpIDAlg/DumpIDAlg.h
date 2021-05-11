#ifndef DumpIDAlg_h
#define DumpIDAlg_h

#include "k4FWCore/DataHandle.h"
#include "GaudiAlg/GaudiAlgorithm.h"
#include "GaudiKernel/NTuple.h"

#include "DetInterface/IGeomSvc.h"

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
    SmartIF<IGeomSvc> m_geosvc;
    dd4hep::Detector* m_dd4hep_geo;
    dd4hep::DDSegmentation::BitFieldCoder* m_decoder;

    DataHandle<edm4hep::SimCalorimeterHitCollection> m_EcalBarrelCol{"EcalBarrelCollection", 
            Gaudi::DataHandle::Reader, this};

private:
    // strore all the id for later analysis
    NTuple::Tuple* m_tuple_id = nullptr;

    NTuple::Item<int> m_id_system;
    NTuple::Item<int> m_id_module;
    NTuple::Item<int> m_id_stave;
    NTuple::Item<int> m_id_tower;
    NTuple::Item<int> m_id_layer;
    NTuple::Item<int> m_id_wafer;
    NTuple::Item<int> m_id_cellX;
    NTuple::Item<int> m_id_cellY;
    
};


#endif
