#ifndef TrackHeedSimTool_h
#define TrackHeedSimTool_h

#include "k4FWCore/DataHandle.h"
#include "GaudiKernel/MsgStream.h"
#include "DetSimInterface/IDedxSimTool.h"
#include <GaudiKernel/AlgTool.h>
#include "edm4hep/MCParticle.h"
#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/SimPrimaryIonizationClusterCollection.h"
#include "TVector3.h"
#include <G4StepPoint.hh>

#include "DD4hep/Segmentations.h"
#include "DD4hep/Printout.h"
#include "DD4hep/Detector.h"
#include "DetInterface/IGeomSvc.h"
#include "DetSegmentation/GridDriftChamber.h"


#include "Garfield/ViewCell.hh"
#include "Garfield/ViewDrift.hh"
#include "Garfield/ViewSignal.hh"
#include "Garfield/ViewMedium.hh"
#include "Garfield/ComponentAnalyticField.hh"
#include "Garfield/MediumMagboltz.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/DriftLineRKF.hh"
#include "Garfield/AvalancheMicroscopic.hh"
#include "Garfield/AvalancheMC.hh"
#include "Garfield/TrackHeed.hh"
#include "Garfield/ComponentNeBem3d.hh"
#include "Garfield/SolidWire.hh"
#include "Garfield/GeometrySimple.hh"
#include "Garfield/MediumConductor.hh"
#include "Garfield/ViewField.hh"

#include <map>
#include <string>

#include "core/session/onnxruntime_cxx_api.h"
#include "core/session/onnxruntime_c_api.h"
using namespace Garfield;

class TrackHeedSimTool: public extends<AlgTool, IDedxSimTool> {
    public:
        using extends::extends;

        StatusCode initialize() override;
        StatusCode finalize() override;
        double dedx(const G4Step* Step) override;
        double dedx(const edm4hep::MCParticle& mc) override;
        double dndx(double betagamma) override;
        void getMom(float ee, float dx, float dy,float dz, float mom[3] );
        void reset(){ 
            m_beginEvt = true;
            m_isFirst = true;
            m_previous_track_ID = 0;
            m_previous_KE = 0;
            m_tot_edep = 0;
            m_tot_length = 0;
        }
        void endOfEvent();
        long long getCellID(float x, float y, float z);
        void wire_xy(float x1, float y1, float z1, float x2, float y2, float z2, float z, float &x, float &y);
        float* NNPred(std::vector<float>& inputs);
        float xy2phi(float x, float y);
        void getLocal(float x1, float y1, float x2, float y2, float& dx, float& dy);
    private:
        //ServiceHandle<IDataProviderSvc> m_eds;
        SmartIF<IGeomSvc> m_geosvc;
        dd4hep::Detector* m_dd4hep; 
        dd4hep::Readout* m_readout;
        dd4hep::DDSegmentation::GridDriftChamber* m_segmentation;
        Gaudi::Property<std::string> m_readout_name{ this, "readout", "DriftChamberHitsCollection"};//readout for getting segmentation
        Gaudi::Property<std::string> m_gas_file{ this, "gas_file", "He_50_isobutane_50.gas"};//gas
        Gaudi::Property<std::string> m_IonMobility{ this, "IonMobility_file", "IonMobility_He+_He.txt"};
        Gaudi::Property<float> m_isob  {this, "isob", 50, ""};
        Gaudi::Property<float> m_he    {this, "he", 50, ""};
        Gaudi::Property<bool> m_debug{this, "debug", false};
        Gaudi::Property<bool> m_use_max_step{this, "use_max_step", false};
        Gaudi::Property<bool> m_update_KE{this, "update_KE", true};
        Gaudi::Property<float> m_max_step   {this, "max_step", 1};//mm
        Gaudi::Property<bool> m_only_primary{this, "only_primary", false};
        Gaudi::Property<bool> m_save_mc{this, "save_mc", false};
        Gaudi::Property<bool> m_save_cellID{this, "save_cellID", true};
        Gaudi::Property<float> m_delta_threshold{this, "delta_threshold", 50};//keV
        Gaudi::Property<float> m_change_threshold {this, "change_threshold", 0.05};
        Gaudi::Property<float> m_BField   {this, "BField", -3};
        Gaudi::Property<float> m_eps     { this, "eps"   , 1e-6  };//very small value, it is returned dedx for unsimulated step (may needed for SimTrackerHit)
        // Output collections
        DataHandle<edm4hep::SimPrimaryIonizationClusterCollection>    m_SimPrimaryIonizationColWriter{"SimPrimaryIonizationClusterCollection", Gaudi::DataHandle::Writer, this};
        edm4hep::SimPrimaryIonizationClusterCollection* m_SimPrimaryIonizationCol;
        // In order to associate MCParticle with contribution, we need to access MC Particle.
        DataHandle<edm4hep::MCParticleCollection> m_mc_handle{"MCParticle", Gaudi::DataHandle::Writer, this};

        TrackHeed* m_track;
        ComponentNeBem3d m_nebem;
        ComponentAnalyticField cmp;
        GeometrySimple m_geo;
        MediumConductor m_metal;
        MediumMagboltz m_gas;
        Sensor* m_sensor;
        std::map<int, std::string> m_particle_map;
        
        int m_previous_track_ID;
        float m_previous_KE;
        int m_current_track_ID;
        int m_current_Parent_ID;
        int m_pdg_code;
        G4StepPoint* m_pre_point;
        G4StepPoint* m_post_point;
        G4double m_total_range;
        bool m_isFirst;
        bool m_beginEvt;
        bool m_change_track;
        edm4hep::MCParticle m_mc_paricle; 
        float m_tot_edep;
        float m_tot_length;
        float m_pa_KE;
  
        G4double m_pre_x  ;
        G4double m_pre_y  ;
        G4double m_pre_z  ;
        G4double m_pre_dx ;
        G4double m_pre_dy ;
        G4double m_pre_dz ;
        G4double m_pre_t  ;
  
        //// sim pulse from NN /// 
        Gaudi::Property<int> m_intra_op_nthreads{ this, "intraOpNumThreads", 1};
        Gaudi::Property<int> m_inter_op_nthreads{ this, "interOpNumThreads", 1};
        std::shared_ptr<Ort::Env> m_env;
        std::shared_ptr<Ort::SessionOptions> m_seesion_options;
        std::shared_ptr<Ort::Session> m_session;
        Ort::AllocatorWithDefaultOptions m_allocator;
        std::vector<const char*> m_input_node_names;
        std::vector<std::vector<int64_t>> m_input_node_dims;
        std::vector<const char*> m_output_node_names;
        std::vector<int64_t> m_output_node_dims;
        #if (ORT_API_VERSION >=13)
        std::vector<Ort::AllocatedStringPtr> m_inputNodeNameAllocatedStrings;
        std::vector<Ort::AllocatedStringPtr> m_outputNodeNameAllocatedStrings;
        #else
        std::vector<const char*> m_inputNodeNameAllocatedStrings;
        std::vector<const char*> m_outputNodeNameAllocatedStrings;
        #endif
  
        Gaudi::Property<bool> m_sim_pulse    { this, "sim_pulse"   , true  };
        Gaudi::Property<std::string> m_model_file{ this, "model", "model_test.onnx"};
        Gaudi::Property<int> m_batchsize     { this, "batchsize", 100};
        Gaudi::Property<float> m_time_scale  { this, "time_scale", 503.0};
        Gaudi::Property<float> m_time_shift  { this, "time_shift", 814.0};
        Gaudi::Property<float> m_amp_scale   { this, "amp_scale" , 1.15 };
        Gaudi::Property<float> m_amp_shift   { this, "amp_shift" , 0.86 };
        Gaudi::Property<float> m_x_scale     { this, "x_scale"   , 9.  };// in mm
        Gaudi::Property<float> m_y_scale     { this, "y_scale"   , 9.  };// in mm
  
  

};

#endif
