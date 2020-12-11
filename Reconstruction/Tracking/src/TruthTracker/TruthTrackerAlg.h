#ifndef TruthTrackerAlg_h
#define TruthTrackerAlg_h

#include "edm4hep/MCParticleCollection.h"
#include "GaudiAlg/GaudiAlgorithm.h"
#include "GaudiKernel/NTuple.h"
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

        virtual StatusCode initialize();
        virtual StatusCode execute();
        virtual StatusCode finalize();
        double cal_dedx_bitrunc(float truncate, std::vector<double> phlist, int & usedhit );
        double BetheBlochEquationDedx(const edm4hep::MCParticle& mcp);

    private:
        SmartIF<IGeomSvc> m_geomSvc;
        dd4hep::Detector* m_dd4hep;
        dd4hep::OverlayedField m_dd4hepField;
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

        Gaudi::Property<int>  m_debug{ this, "debug", false};
        Gaudi::Property<float>  m_truncate{ this, "truncate", 0.7};
        Gaudi::Property<float>  m_mom_resolution{ this, "mom_resolution", 0.};
        Gaudi::Property<bool>  m_WriteAna{ this, "WriteAna", false};

        Gaudi::Property<float> m_material_Z{this, "material_Z", 7};//Default is Nitrogen
        Gaudi::Property<float> m_material_A{this, "material_A", 14};
        Gaudi::Property<float> m_track_dedx_resolution{this, "track_dedx_resolution", 0.};
        Gaudi::Property<float> m_dedx_scale{this, "dedx_scale", 1};
        Gaudi::Property<float> m_dedx_resolution{this, "dedx_resolution", 0.};
        float m_me;// Here me is the electron rest mass
        float m_K; // K was set as a constant.

        NTuple::Tuple* m_tuple = nullptr ;
        NTuple::Item<long>   m_hit;
        NTuple::Item<double>   m_track_dedx;
        NTuple::Item<double>   m_track_dedx_BB;
        NTuple::Item<double>   m_track_px;
        NTuple::Item<double>   m_track_py;
        NTuple::Item<double>   m_track_pz;
        NTuple::Item<double>   m_track_mass;
        NTuple::Item<int>   m_track_pid;
        NTuple::Array<double  > m_hit_x   ;
        NTuple::Array<double  > m_hit_y   ;
        NTuple::Array<double  > m_hit_z   ;
        NTuple::Array<double  > m_hit_dedx;
        NTuple::Item<double>   m_mc_px;
        NTuple::Item<double>   m_mc_py;
        NTuple::Item<double>   m_mc_pz;

};

#endif
