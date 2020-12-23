#ifndef G2CDARBORALG_H
#define G2CDARBORALG_H

#include "k4FWCore/DataHandle.h"
#include "GaudiAlg/GaudiAlgorithm.h"
#include "Gaudi/Property.h"
#include "edm4hep/EventHeader.h"
#include "edm4hep/EventHeaderCollection.h"
#include "edm4hep/SimCalorimeterHitConst.h"
#include "edm4hep/SimCalorimeterHit.h"
#include "edm4hep/CalorimeterHit.h"
#include "edm4hep/CalorimeterHitCollection.h"
#include "edm4hep/SimCalorimeterHitCollection.h"
#include "edm4hep/MCRecoCaloAssociationCollection.h"
#include "edm4hep/MCParticleCollection.h"

#include <DDRec/DetectorData.h>
#include <DDRec/CellIDPositionConverter.h>
#include "DetInterface/IGeomSvc.h"

#include <string>
#include <iostream>
#include <fstream>

#include <TObject.h>
#include <TNtuple.h>
#include <TTree.h>
#include <TFile.h>
#include <TF1.h>
#include <TH1F.h>
class TTree;

class G2CDArborAlg  : public GaudiAlgorithm
{
public:
     /* Processor*  newProcessor() { return new G2CDArborAlg ; } */

     G2CDArborAlg(const std::string& name, ISvcLocator* svcLoc);
     /* G2CDArborAlg(); */
     /* ~G2CDArborAlg() {}; */

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
     std::string GetLayerCoding(const std::string &encodingString) const;

     /* DataHandle<edm4hep::EventHeaderCollection> m_headerCol{"EventHeaderCol", Gaudi::DataHandle::Reader, this}; */

     // Input collections
     DataHandle<edm4hep::MCParticleCollection> r_mcParticle{"MCParticle", Gaudi::DataHandle::Reader, this};

     typedef DataHandle<edm4hep::SimCalorimeterHitCollection>  SimCaloType;

     Gaudi::Property<std::vector<std::string>> m_ecalColNames{this, "ECALCollections", {"EcalBarrelSiliconCollection", "EcalEndcapSiliconCollection", "EcalEndcapRingCollection"}, "Input ECAL Hits Collection Names"};
     std::vector<SimCaloType*> _ecalCollections;

     Gaudi::Property<std::vector<std::string>> m_ecalPreShowerColNames{this, "HcalHitCollections", {"EcalBarrelSiliconPreShowerCollection", "EcalEndcapRingPreShowerCollection", "EcalEndcapSiliconPreShowerCollection"}, "Hit Collection Names"};
     std::vector<SimCaloType*> _EcalPreShowerCollections;

     Gaudi::Property<std::vector<std::string>> m_hcalColNames{this, "HCALCollections", {"HcalBarrelRegCollection", "HcalEndcapsCollection", "HcalEndcapRingsCollection"}, "HCAL Collection Names"};
     std::vector<SimCaloType*> _hcalCollections;

     Gaudi::Property<std::vector<std::string>> m_ecalReadoutNames{this, "ECALReadOutNames", {"EcalBarrelCollection", "EcalEndcapsCollection","EcalEndcapRingCollection"}, "Name of readouts"};
     Gaudi::Property<std::vector<std::string>> m_hcalReadoutNames{this, "HCALReadOutNames", {"HcalBarrelCollection", "HcalEndcapsCollection","HcalEndcapRingCollection"}, "Name of readouts"};
     std::map<std::string, std::string> m_col_readout_map;

     // Output collections
     typedef DataHandle<edm4hep::CalorimeterHitCollection>  CaloHitType;

     Gaudi::Property<std::vector<std::string>> m_ecalOutputColNames{this, "DigiECALCollection", {"ECALBarrel", "ECALEndcap", "ECALOther"}, "Name of Digitized ECAL Hit Collections"};
     std::vector<CaloHitType*> _outputEcalCollections;

     Gaudi::Property<std::vector<std::string>> m_hcalOutputColNames{this, "DigiHCALCollection", {"HCALBarrel", "HCALEndcap", "HCALOther"}, "Name of Digitized HCAL Hit Collections"};
     std::vector<CaloHitType*> _outputHcalCollections;

     /* typedef DataHandle<edm4hep::MCRecoCaloAssociationCollection>  McRecoCaloAssoType; */
     /* Gaudi::Property<std::vector<std::string>> m_caloTruthLinkColName{this, "caloTruthLinkCollection", {}, "caloTruthLinkCollection"}; */
     /* std::vector<McRecoCaloAssoType*> _caloTruthLinkCollection; */
     DataHandle<edm4hep::MCRecoCaloAssociationCollection> _caloTruthLinkCollection{"MCRecoCaloAssociationCollection", Gaudi::DataHandle::Writer, this};

     mutable Gaudi::Property<std::vector<float>> m_ChargeSpatialDistri{this, "ChargeSpatialDistribution", {0.1, 0.2, 0.4, 0.2, 0.1}, "Spactial Distribution of MIP charge X*Y;"};
     mutable Gaudi::Property<std::vector<float>> _calibCoeffEcal{this, "CalibrECAL", {40.91, 81.81}, "Calibration coefficients for ECAL"};
     mutable Gaudi::Property<std::vector<float>> _ShowerPositionShiftID{this, "PositionShiftID", {0, 0, 0}, "Global Position Shift For Overlay"}; //should be of the form deltaI, J, K */

     Gaudi::Property<bool>   m_readLCIO{this, "ReadLCIO", false, "Read sim file with LCIO"};
     Gaudi::Property<int>    m_reportEvery{this, "EventReportEvery", 100, "Event ID report every"};

     Gaudi::Property<int>    _NEcalThinLayer{this, "NumThinEcalLayer", 20, "Num of thiner Ecal layers"};
     Gaudi::Property<float>  _thresholdEcal{this, "ECALThreshold", (float)5.0e-5, "Threshold for ECAL Hits in GeV"};
     Gaudi::Property<float>  _thresholdHcal{this, "HCALThreshold", (float)0.11, "Threshold for HCAL Hits in GeV"};
     Gaudi::Property<int>    _DigiCellSize{this, "DigiCellSize", 10, "Size of Digitized Cell (in mm)"};
     Gaudi::Property<float>  _ShiftInX{this, "ShiftInX", (float)0.0, "Shift Distance in X directoin (in mm) NP only"};
     Gaudi::Property<int>    _UsingDefaultDetector{this, "UsingDefaultDetector", 0, "Flag Parameter Setting (0 ~ self definition, 1 ~ MircoMegas, 2 ~ GRPC_PS, 3 ~ GRPC_SPS)"};
     Gaudi::Property<float>  _PolyaParaA{this, "PolyaParaA", 0.7, "Polya: x^A*exp(-b*x) + c"};
     Gaudi::Property<float>  _PolyaParaB{this, "PolyaParaB", 0.045, "Polya: x^a*exp(-B*x) + c"};
     Gaudi::Property<float>  _PolyaParaC{this, "PolyaParaC", 0.03, "Polya: x^a*exp(-b*x) + C"};
     Gaudi::Property<float>  _ChanceOfKink{this, "ChanceOfKink", 0.0, "Chance of one boundary hit to create a multiple hit with boosted charge"};
     Gaudi::Property<float>  _KinkHitChargeBoost{this, "KinkHitChargeBoost", 1.0, "Scale factor of Charge on boosted multiple hits"};

     /* std::string _treeFileName; */
     /* std::string _treeName; */
     /* std::string _colName; */
     /* std::vector<std::string> _caloTruthLinkCollection; */
     /* std::vector<std::string> _hcalCollections; */
     /* std::vector<std::string> _outputHcalCollections; */
     /* std::vector<std::string> _ecalCollections; */
     /* std::vector<std::string> _outputEcalCollections; */
     /* std::vector<std::string> _EcalPreShowerCollections; */
     std::vector<float> _ChargeSpatialDistri;
     /* std::vector<float> _calibCoeffEcal; */
     /* std::vector<float> _ShowerPositionShiftID; //should be of the form deltaI, J, K */

     std::map <int, std::pair<float, float> >WeightVector;

     const double m_pi;

     /* float _thresholdEcal; */
     /* float _thresholdHcal; */
     /* int _NEcalThinLayer;  */
     int _overwrite;
     /* int _DigiCellSize;  */
     /* int _UsingDefaultDetector;  */
     /* float _PolyaParaA, _PolyaParaB, _PolyaParaC;  */
     /* float _ChanceOfKink, _KinkHitChargeBoost;  */
     TTree *_outputTree;
     TH1F *_NH1stLayer, *_NH8thLayer; 
     TF1 * _QPolya; 

     int _Num;
     int _eventNr; 

     int _M, _S, _I, _J, _K, _Seg;
     float _PosX, _PosY, _PosZ; 
     /* float _EDepo, _Charge, _ShiftInX;  */
     int _NHit1mm, _NHit1mmCenter, _NHit1mmCorner, _NHit1mmSide;
     int _TotalNHit1mm, _TotalNHit, _TotalMultiHit; 
     int _N1, _N2, _N3; 

     std::string _fileName;
     std::ostream *_output;
     std::string m_encoder_str;

     SmartIF<IGeomSvc> m_geosvc;
     dd4hep::Detector* m_dd4hep_geo;
     dd4hep::DDSegmentation::BitFieldCoder* m_decoder;
};

#endif


