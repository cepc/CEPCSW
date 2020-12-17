#ifndef RecDCHDedxAlg_h
#define RecDCHDedxAlg_h

#include "edm4hep/MCParticleCollection.h"
#include "GaudiAlg/GaudiAlgorithm.h"
#include "GaudiKernel/NTuple.h"
#include "k4FWCore/DataHandle.h"
#include "DD4hep/Fields.h"
#include "DetSimInterface/IDedxSimTool.h"
#include "GaudiKernel/ToolHandle.h"

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

class RecDCHDedxAlg: public GaudiAlgorithm
{
    public:
        RecDCHDedxAlg(const std::string& name, ISvcLocator* svcLoc);

        virtual StatusCode initialize();
        virtual StatusCode execute();
        virtual StatusCode finalize();
        double cal_dedx_bitrunc(float truncate, std::vector<double> phlist, int & usedhit );

    private:
        ToolHandle<IDedxSimTool> m_dedx_simtool;
        Gaudi::Property<std::string>  m_dedx_sim_option{ this, "sampling_option", ""};

        //reader
        DataHandle<edm4hep::MCRecoTrackerAssociationCollection> m_dcHitAssociationCol{ "DCHitAssociationCollection",Gaudi::DataHandle::Reader, this};
        //writer
        DataHandle<edm4hep::TrackCollection> m_dcTrackCol{"DCTrackCollection", Gaudi::DataHandle::Writer, this};

        Gaudi::Property<int>  m_debug{ this, "debug", false};
        Gaudi::Property<int>  m_method{ this, "method", 1};
        Gaudi::Property<float>  m_truncate{ this, "truncate", 0.7};
        Gaudi::Property<bool>  m_WriteAna{ this, "WriteAna", false};
        Gaudi::Property<float> m_scale{this, "dedx_scale", 1};
        Gaudi::Property<float> m_resolution{this, "dedx_resolution", 0.};

        NTuple::Tuple* m_tuple = nullptr ;
        NTuple::Item<long>   m_n_track;
        NTuple::Item<long>   m_hit;
        NTuple::Array<double>   m_track_dedx;
        NTuple::Array<double>   m_track_dedx_BB;
        NTuple::Array<double>   m_track_px;
        NTuple::Array<double>   m_track_py;
        NTuple::Array<double>   m_track_pz;
        NTuple::Array<double>   m_track_mass;
        NTuple::Array<int>      m_track_pid;
        NTuple::Array<double> m_hit_x   ;
        NTuple::Array<double> m_hit_y   ;
        NTuple::Array<double> m_hit_z   ;
        NTuple::Array<double> m_hit_dedx;
        NTuple::Item<double>  m_mc_px;
        NTuple::Item<double>  m_mc_py;
        NTuple::Item<double>  m_mc_pz;

};

#endif
