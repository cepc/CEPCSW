#ifndef DCH_DIGI_ALG_H
#define DCH_DIGI_ALG_H

#include "FWCore/DataHandle.h"
#include "GaudiAlg/GaudiAlgorithm.h"
#include "edm4hep/SimTrackerHitCollection.h"
#include "edm4hep/TrackerHitCollection.h"
#include "edm4hep/MCRecoTrackerAssociationCollection.h"

#include <DDRec/DetectorData.h>
#include <DDRec/CellIDPositionConverter.h>
#include "DetInterface/IGeomSvc.h"
#include "DDSegmentation/Segmentation.h"
#include "DetSegmentation/GridDriftChamber.h"

#include "TVector3.h"
#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"



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
  void Reset();
 
protected:

  SmartIF<IGeomSvc> m_geosvc;
  typedef std::vector<float> FloatVec;
  int _nEvt ;

  TFile* m_fout;
  TTree* m_tree;
  std::vector<int  > m_chamber   ;
  std::vector<int  > m_layer     ;
  std::vector<int  > m_cell      ;
  std::vector<float> m_cell_x    ;
  std::vector<float> m_cell_y    ;
  std::vector<float> m_simhit_x  ;
  std::vector<float> m_simhit_y  ;
  std::vector<float> m_simhit_z  ;
  std::vector<float> m_hit_x     ;
  std::vector<float> m_hit_y     ;
  std::vector<float> m_hit_z     ;
  std::vector<float> m_dca       ;
  std::vector<float> m_hit_dE    ;
  std::vector<float> m_hit_dE_dx ;



  dd4hep::rec::CellIDPositionConverter* m_cellIDConverter;
  dd4hep::DDSegmentation::GridDriftChamber* m_segmentation;
  dd4hep::DDSegmentation::BitFieldCoder* m_decoder;
  
  Gaudi::Property<std::string> m_readout_name{ this, "readout", "DriftChamberHitsCollection"};//readout for getting segmentation
 
  Gaudi::Property<float> m_res_x     { this, "res_x", 0.11};//mm
  Gaudi::Property<float> m_res_y     { this, "res_y", 0.11};//mm
  Gaudi::Property<float> m_res_z     { this, "res_z", 1   };//mm
  Gaudi::Property<float> m_velocity  { this, "drift_velocity", 40};// um/ns
  Gaudi::Property<float> m_mom_threshold { this, "mom_threshold", 0};// GeV
  Gaudi::Property<bool>  m_WriteAna { this, "WriteAna", false};

  Gaudi::Property<std::string> m_Output { this, "output", "ana_DCH_digi.root"};

  // Input collections
  DataHandle<edm4hep::SimTrackerHitCollection> r_SimDCHCol{"DriftChamberHitsCollection", Gaudi::DataHandle::Reader, this};
  // Output collections
  DataHandle<edm4hep::TrackerHitCollection>    w_DigiDCHCol{"DigiDCHitsCollection", Gaudi::DataHandle::Writer, this};
  DataHandle<edm4hep::MCRecoTrackerAssociationCollection>    w_AssociationCol{"DCHitAssociationCollection", Gaudi::DataHandle::Writer, this};
};

#endif
