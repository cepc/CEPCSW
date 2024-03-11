#ifndef DCH_DIGI_ALG_H
#define DCH_DIGI_ALG_H

#include "k4FWCore/DataHandle.h"
#include "GaudiKernel/NTuple.h"
#include "GaudiAlg/GaudiAlgorithm.h"
#include "edm4hep/SimTrackerHitCollection.h"
#include "edm4hep/TrackerHitCollection.h"
#include "edm4hep/MCRecoTrackerAssociationCollection.h"
#include "edm4hep/MCParticle.h"

#include <DDRec/DetectorData.h>
#include <DDRec/CellIDPositionConverter.h>
#include "DetInterface/IGeomSvc.h"
#include "DDSegmentation/Segmentation.h"
#include "DetSegmentation/GridDriftChamber.h"

#include "TVector3.h"
#include "TRandom3.h"

class DCHDigiAlg : public GaudiAlgorithm
{
 
public:
 
  DCHDigiAlg(const std::string& name, ISvcLocator* svcLoc);
 
  /** Called at the begin of the job before anything is read.
   * Use to initialize the processor, e.g. book histograms.
   */
  virtual StatusCode initialize() ;
 
  /** Called for every event - the working horse.
   */
  virtual StatusCode execute() ; 
 
  /** Called after data processing for clean up.
   */
  virtual StatusCode finalize() ;
 
protected:

  SmartIF<IGeomSvc> m_geosvc;
  typedef std::vector<float> FloatVec;
  int _nEvt ;

  TRandom3 fRandom;

  NTuple::Tuple* m_tuple = nullptr ;
  NTuple::Item<int> m_evt;
  NTuple::Item<long>   m_n_sim;
  NTuple::Item<long>   m_n_digi;
  NTuple::Item<float> m_time;
  NTuple::Array<int  > m_chamber   ;
  NTuple::Array<int  > m_layer     ;
  NTuple::Array<int  > m_cell      ;
  NTuple::Array<float> m_cell_x    ;
  NTuple::Array<float> m_cell_y    ;
  NTuple::Array<float> m_cell1_x   ;
  NTuple::Array<float> m_cell1_y   ;
  NTuple::Array<float> m_simhit_x  ;
  NTuple::Array<float> m_simhit_y  ;
  NTuple::Array<float> m_simhit_z  ;
  NTuple::Array<float> m_simhitT ;
  NTuple::Array<float> m_simhitmom  ;
  NTuple::Array<int> m_simPDG  ;
  NTuple::Array<float> m_hit_x     ;
  NTuple::Array<float> m_hit_y     ;
  NTuple::Array<float> m_hit_z     ;
  NTuple::Array<float> m_mom_x     ;
  NTuple::Array<float> m_mom_y     ;
  NTuple::Array<float> m_dca       ;
  NTuple::Array<float> m_Simdca       ;
  NTuple::Array<float> m_poca_x    ;
  NTuple::Array<float> m_poca_y    ;
  NTuple::Array<float> m_hit_dE    ;
  NTuple::Array<float> m_hit_dE_dx ;
  NTuple::Array<double> m_truthlength ;

  clock_t m_start,m_end;

  dd4hep::rec::CellIDPositionConverter* m_cellIDConverter;
  dd4hep::DDSegmentation::GridDriftChamber* m_segmentation;
  dd4hep::DDSegmentation::BitFieldCoder* m_decoder;
  
  Gaudi::Property<std::string> m_readout_name{ this, "readout", "DriftChamberHitsCollection"};//readout for getting segmentation
 
  Gaudi::Property<float> m_res_x     { this, "res_x", 0.11};//mm
  Gaudi::Property<float> m_res_y     { this, "res_y", 0.11};//mm
  Gaudi::Property<float> m_res_z     { this, "res_z", 1   };//mm
  Gaudi::Property<float> m_velocity  { this, "drift_velocity", 40};// um/ns
  Gaudi::Property<float> m_mom_threshold { this, "mom_threshold", 0};// GeV
  Gaudi::Property<float> m_mom_threshold_high { this, "mom_threshold_high", 1e9};// GeV
  Gaudi::Property<float> m_edep_threshold{ this, "edep_threshold", 0};// GeV
  Gaudi::Property<bool>  m_WriteAna { this, "WriteAna", false};
  Gaudi::Property<bool>  m_debug{ this, "debug", false};
  Gaudi::Property<double>  m_wireEff{ this, "wireEff", 1.0};

  // Input collections
  DataHandle<edm4hep::SimTrackerHitCollection> r_SimDCHCol{"DriftChamberHitsCollection", Gaudi::DataHandle::Reader, this};
  // Output collections
  DataHandle<edm4hep::TrackerHitCollection>    w_DigiDCHCol{"DigiDCHitCollection", Gaudi::DataHandle::Writer, this};
  DataHandle<edm4hep::MCRecoTrackerAssociationCollection>    w_AssociationCol{"DCHitAssociationCollection", Gaudi::DataHandle::Writer, this};
};

#endif
