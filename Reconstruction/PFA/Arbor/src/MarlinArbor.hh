#ifndef _MarlinArbor_hh_
#define _MarlinArbor_hh_

#include <string>
#include <iostream>
#include <fstream>
#include "k4FWCore/DataHandle.h"
#include "GaudiAlg/GaudiAlgorithm.h"
#include "Gaudi/Property.h"
#include "edm4hep/EventHeader.h"
#include "edm4hep/EventHeaderCollection.h"
#include "edm4hep/SimCalorimeterHitConst.h"
#include "edm4hep/SimCalorimeterHit.h"
#include "edm4hep/CalorimeterHit.h"
#include "edm4hep/CalorimeterHitCollection.h"
#include "edm4hep/Cluster.h"
#include "edm4hep/ClusterCollection.h"
#include "edm4hep/SimCalorimeterHitCollection.h"
#include "edm4hep/MCRecoCaloAssociationCollection.h"
#include "edm4hep/MCParticleCollection.h"

#include <DDRec/DetectorData.h>
#include <DDRec/CellIDPositionConverter.h>
#include "DetInterface/IGeomSvc.h"

#include "Arbor.h"
#include "ArborToolLCIO.hh"

#include <TTree.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>

class TTree;

//class CollectionMaps
//{
//public:
//    CollectionMaps();
//    void clear();
//    std::map<std::string, std::vector<edm4hep::MCParticle> >     collectionMap_MC;
//    std::map<std::string, std::vector<edm4hep::CalorimeterHit> > collectionMap_CaloHit;
//};

class MarlinArbor  : public GaudiAlgorithm
{
	public:

		MarlinArbor(const std::string& name, ISvcLocator* svcLoc);


		virtual StatusCode initialize() ;

		virtual StatusCode execute() ;

		virtual StatusCode finalize() ;
		void HitsPreparation(); // LCEvent * evtP);

		void MakeIsoHits(std::vector<edm4hep::ConstCalorimeterHit> inputCaloHits, DataHandle<edm4hep::CalorimeterHitCollection>& m_isohitcol);
	protected:
		std::string _colName;
		std::vector<std::string> _CalCollections;
		std::vector<std::string> _SimCalCollections;
		std::vector<std::string> _garlicCollections;
		std::vector<std::string> _endHitTrackCollections;
		std::vector<std::string> _EcalPreShowerCollections;
		std::vector<std::string> _EcalCalCollections;
		std::vector<std::string> _HcalCalCollections;
		std::vector<float> _cepc_thresholds; 


	     typedef DataHandle<edm4hep::CalorimeterHitCollection>  CaloType;
	     Gaudi::Property<std::vector<std::string>> m_ecalColNames{this, "ECALCollections", {"ECALBarrel", "ECALEndcap", "ECALOther"}, "Input ECAL Hits Collection Names"};
	     Gaudi::Property<std::vector<std::string>> m_hcalColNames{this, "HCALCollections", {"HCALBarrel", "HCALEndcap", "HCALOther"}, "Input HCAL Hits Collection Names"};


     Gaudi::Property<std::vector<std::string>> m_ecalReadoutNames{this, "ECALReadOutNames", {"EcalBarrelCollection", "EcalEndcapsCollection","EcalEndcapRingCollection"}, "Name of readouts"};
     Gaudi::Property<std::vector<std::string>> m_hcalReadoutNames{this, "HCALReadOutNames", {"HcalBarrelCollection", "HcalEndcapsCollection","HcalEndcapRingCollection"}, "Name of readouts"};
     std::map<std::string, std::string> m_col_readout_map;

	     std::vector<CaloType*> _ecalCollections;
	     std::vector<CaloType*> _hcalCollections;
	     std::vector<CaloType*> _calCollections;

		TTree *_outputTree;
		std::string _treeFileName; 
		int _EH; 
		float _HitPos[3];
		float _BushP[3];
		float _CloseDis; 
		float _HitEnergy;

		int _CellSize; 
		int _CaloTrackLengthCut;

		int _Num, _Seg, _eventNr;
		int numElements;

		float _DHCALFirstThreshold; 
		float _InitLinkDisThreshold;

		bool _DHCALSimuDigiMode; 
		bool _FlagInputSimHit;
		bool _FlagMutePhoton;
		bool _FlagMuteChargeParticle;
		bool _FlagMuteGarlicHits;
		bool _FlagUseTrackerEndHit; 
		std::string m_encoder_str;

		DataHandle<edm4hep::ClusterCollection>   branchCol{"EHBushes",Gaudi::DataHandle::Writer, this};
		DataHandle<edm4hep::CalorimeterHitCollection>   m_isohitcol{"IsoHits",Gaudi::DataHandle::Writer, this};
		TH2F *_h1, *_h2, *_h7; 
		TH1F *_h3, *_h4, *_h5, *_h6; 
		std::ostream _output;
		float _HLayerCut;

     SmartIF<IGeomSvc> m_geosvc;
     dd4hep::DDSegmentation::BitFieldCoder* m_decoder;
     ArborToolLCIO * m_ArborToolLCIO;
};


#endif


