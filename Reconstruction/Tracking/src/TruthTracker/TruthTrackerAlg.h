#ifndef TruthTrackerAlg_h
#define TruthTrackerAlg_h

#include "GaudiAlg/GaudiAlgorithm.h"
#include "k4FWCore/DataHandle.h"
#include "DD4hep/Fields.h"

class IGeomSvc;
namespace dd4hep {
    class Detector;
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

        virtual StatusCode initialize() override;
        virtual StatusCode execute() override;
        virtual StatusCode finalize() override;

    private:
        SmartIF<IGeomSvc> m_geomSvc;
        dd4hep::Detector* m_dd4hep;
        dd4hep::OverlayedField m_dd4hepField;
        dd4hep::DDSegmentation::GridDriftChamber* m_gridDriftChamber;
        dd4hep::DDSegmentation::BitFieldCoder* m_decoder;

        //reader
        DataHandle<edm4hep::MCParticleCollection> m_mcParticleCol{
            "MCParticle", Gaudi::DataHandle::Reader, this};
        DataHandle<edm4hep::TrackerHitCollection> m_DCDigiCol{
            "DigiDCHitCollection", Gaudi::DataHandle::Reader, this};
        DataHandle<edm4hep::MCRecoTrackerAssociationCollection>
            m_DCHitAssociationCol{ "DCHitAssociationCollection",
                Gaudi::DataHandle::Reader, this};
        DataHandle<edm4hep::TrackCollection>
            m_siSubsetTrackCol{ "SiSubsetTrackCollection",
                Gaudi::DataHandle::Reader, this};
        //writer
        DataHandle<edm4hep::TrackCollection> m_DCTrackCol{
            "DCTrackCollection", Gaudi::DataHandle::Writer, this};
        DataHandle<edm4hep::TrackCollection> m_SDTTrackCol{
            "SDTTrackCollection", Gaudi::DataHandle::Writer, this};
        DataHandle<edm4hep::ReconstructedParticleCollection> m_DCRecParticleCol{
            "DCRecParticleCollection", Gaudi::DataHandle::Writer, this};
        DataHandle<edm4hep::MCRecoParticleAssociationCollection>
            m_DCRecParticleAssociationCol{
                "DCRecMCRecoParticleAssociationCollection",
                Gaudi::DataHandle::Writer, this};

        //readout for getting segmentation
        Gaudi::Property<std::string> m_readout_name{this, "readout",
            "DriftChamberHitsCollection"};
        Gaudi::Property<bool> m_writeRecParticle{this,"writeRecParticle",false};
        Gaudi::Property<float> m_resPT{this,"resPT",0};//ratio
        Gaudi::Property<float> m_resPz{this,"resPz",0};//ratio
        Gaudi::Property<float> m_resMomPhi{this,"resMomPhi",0};//radian
        Gaudi::Property<float> m_resMomTheta{this,"resMomTheta",0};//radian
        Gaudi::Property<float> m_resVertexX{this,"resVertexX",0.003};//3um
        Gaudi::Property<float> m_resVertexY{this,"resVertexY",0.003};//3um
        Gaudi::Property<float> m_resVertexZ{this,"resVertexZ",0.003};//3um
        Gaudi::Property<int> m_maxDCDigiCut{this,"maxDigiCut",1e6};
};

#endif
