#ifndef CRD_HCAL_DIGI_ALG_H 
#define CRD_HCAL_DIGI_ALG_H

#include "k4FWCore/DataHandle.h"
#include "GaudiAlg/GaudiAlgorithm.h"
#include "edm4hep/MutableCaloHitContribution.h"
#include "edm4hep/MutableSimCalorimeterHit.h"
#include "edm4hep/SimCalorimeterHit.h"
#include "edm4hep/CalorimeterHit.h"
#include "edm4hep/CalorimeterHitCollection.h"
#include "edm4hep/SimCalorimeterHitCollection.h"
#include "edm4hep/MCRecoCaloAssociationCollection.h"
#include "edm4hep/MCRecoCaloParticleAssociationCollection.h"
#include "edm4hep/MCParticle.h"

#include <DDRec/DetectorData.h>
#include <DDRec/CellIDPositionConverter.h>
#include <DD4hep/Segmentations.h> 
#include "DetInterface/IGeomSvc.h"
#include "TVector3.h"
#include "TRandom3.h"
#include "TFile.h"
#include "TString.h"
#include "TH3.h"
#include "TH1.h"

#include <cstdlib>
#include "time.h"
#include <TTimeStamp.h> 
#include <ctime>

#define C 299.79  // unit: mm/ns
#define PI 3.141592653

class CRDHcalDigiAlg : public GaudiAlgorithm
{
 
public:
 
  CRDHcalDigiAlg(const std::string& name, ISvcLocator* svcLoc);
 
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

	void Clear();

protected:

  SmartIF<IGeomSvc> m_geosvc;
  //SmartIF<ICRDEcalSvc> m_edmsvc;
  typedef std::vector<float> FloatVec;
  typedef std::map<const edm4hep::MCParticle, float> MCParticleToEnergyWeightMap;

	int _nEvt ;
	TRandom3 rndm;
	TFile* m_wfile;
	TTree* t_simHit;

	FloatVec m_simHit_x, m_simHit_y, m_simHit_z, m_simHit_E, m_simHit_slice, m_simHit_layer, m_simHit_tower, m_simHit_stave, m_simHit_module;
  std::vector<unsigned long long> m_simHit_cellID;
  std::vector<int> m_simHit_steps;


	dd4hep::rec::CellIDPositionConverter* m_cellIDConverter;
	dd4hep::DDSegmentation::BitFieldCoder* m_decoder;

	Gaudi::Property<float> m_scale{ this, "Scale", 1 };

  // Input collections
  DataHandle<edm4hep::SimCalorimeterHitCollection> r_SimCaloCol{"SimCaloCol", Gaudi::DataHandle::Reader, this};
  mutable Gaudi::Property<std::string> _readoutName{this, "ReadOutName", "CaloHitsCollection", "Readout name"};
  mutable Gaudi::Property<std::string> _filename{this, "OutFileName", "testout.root", "Output file name"};

  //Input parameters
  mutable Gaudi::Property<int>   _writeNtuple{this,  "WriteNtuple", 1, "Write ntuple"};
  mutable Gaudi::Property<int>   _Nskip{this,  "SkipEvt", 0, "Skip event"};
  mutable Gaudi::Property<float> _seed{this,   "Seed", 2131, "Random Seed"};
  mutable Gaudi::Property<int>  _Debug{this,   "Debug", 0, "Debug level"};
  mutable Gaudi::Property<float> _MIPCali {this,   "MIPResponse",  0.0005, "MIP response (/GeV)"};
  mutable Gaudi::Property<float> _Eth_Mip {this,   "MIPThreshold", 0.5, "Energy Threshold (/MIP)"};
  mutable Gaudi::Property<float> r_cali{this,  "CalibrHCAL", 1, "Global calibration coefficients for HCAL"};


  // Output collections
  DataHandle<edm4hep::CalorimeterHitCollection>    w_DigiCaloCol{"DigiCaloCol", Gaudi::DataHandle::Writer, this};
  DataHandle<edm4hep::MCRecoCaloAssociationCollection>    w_CaloAssociationCol{"HCALBarrelAssoCol", Gaudi::DataHandle::Writer, this};
  DataHandle<edm4hep::MCRecoCaloParticleAssociationCollection>    w_MCPCaloAssociationCol{"HCALBarrelParticleAssoCol", Gaudi::DataHandle::Writer, this};
};

#endif
