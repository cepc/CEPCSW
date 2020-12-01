#ifndef TruthTrackerAlg_h
#define TruthTrackerAlg_h

#include "GaudiAlg/GaudiAlgorithm.h"
#include "k4FWCore/DataHandle.h"
#include "GaudiKernel/NTuple.h"
#include "DD4hep/Fields.h"

class IGeomSvc;
namespace dd4hep {
    class Detector;
    //class rec::CellIDPositionConverter;
    namespace DDSegmentation{
        class GridDriftChamber;
        class BitFieldCoder;
    }
}
namespace edm4hep {
    class MCParticleCollection;
    class TrackerHitCollection;
    class TrackCollection;
    class MCRecoTrackerAssociationCollection;
    class ReconstructedParticleCollection;
    class MCRecoParticleAssociationCollection;
}

class TruthTrackerAlg: public GaudiAlgorithm
{
    public:
        TruthTrackerAlg(const std::string& name, ISvcLocator* svcLoc);

        virtual StatusCode initialize();
        virtual StatusCode execute();
        virtual StatusCode finalize();

    private:
        SmartIF<IGeomSvc> m_geomSvc;
        dd4hep::Detector* m_dd4hep;
        dd4hep::OverlayedField m_dd4hepField;
        //dd4hep::rec::CellIDPositionConverter* m_cellIDConverter;
        dd4hep::DDSegmentation::GridDriftChamber* m_gridDriftChamber;
        dd4hep::DDSegmentation::BitFieldCoder* m_decoder;

        //reader
        DataHandle<edm4hep::MCParticleCollection> m_mcParticleCol{
            "MCParticle", Gaudi::DataHandle::Reader, this};
        DataHandle<edm4hep::TrackerHitCollection> m_digiDCHitsCol{
            "DigiDCHitsCollection", Gaudi::DataHandle::Reader, this};
        DataHandle<edm4hep::MCRecoTrackerAssociationCollection>
            m_dcHitAssociationCol{ "DCHitAssociationCollection",
                Gaudi::DataHandle::Reader, this};
        //writer
        DataHandle<edm4hep::TrackCollection> m_dcTrackCol{"DCTrackCollection",
            Gaudi::DataHandle::Writer, this};
        DataHandle<edm4hep::ReconstructedParticleCollection> m_dcRecParticleCol{
            "DCRecParticleCollection", Gaudi::DataHandle::Writer, this};
        DataHandle<edm4hep::MCRecoParticleAssociationCollection>
            m_dcRecParticleAssociationCol{"DCRecMCRecoParticleAssociationCollection",
                Gaudi::DataHandle::Writer, this};

        //readout for getting segmentation
        Gaudi::Property<std::string> m_readout_name{this, "readout",
            "DriftChamberHitsCollection"};

        int m_debug;
        //strore all the id for later analysis
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
