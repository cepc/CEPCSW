#ifndef _BushConnect_hh_
#define _BushConnect_hh_

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


class BushConnect  : public GaudiAlgorithm
{
	public:

		BushConnect(const std::string& name, ISvcLocator* svcLoc);


		virtual StatusCode initialize() ;

		virtual StatusCode execute() ;

		virtual StatusCode finalize() ;
		void Clean(); 
		void TrackSort(); 
		void BushSelfMerge(); 
		void TagCore(); 
		void ParticleReco(); 


	protected:

		std::vector<edm4hep::ConstCluster> SortedSMBushes;
		std::vector<edm4hep::ConstTrack> SortedTracks;
		std::map<edm4hep::ConstTrack, float> Track_Energy;
		std::map<edm4hep::ConstTrack, TVector3> Track_P3;
		std::map<edm4hep::ConstTrack, int> Track_Type;
		std::map<edm4hep::ConstTrack, float> Track_Theta;
		std::map<edm4hep::ConstTrack, float> Track_Phi;	

		std::map<edm4hep::ConstCluster, int> ClusterType_1stID;
		std::map<edm4hep::ReconstructedParticle, int> ChCoreID; 

		std::vector<edm4hep::ConstCluster> ecalchcore_tight;         //TightCores
		std::vector<edm4hep::ConstCluster> ecalchcore_medium;
		std::vector<edm4hep::ConstCluster> ecalchcore_loose;         //LooseCores    Let's also get
		std::vector<edm4hep::ConstCluster> ecalchcore; 		    //Above three
		std::vector<edm4hep::ConstCluster> ecalnecore;
		std::vector<edm4hep::ConstCluster> ecalnecore_EM;
		std::vector<edm4hep::ConstCluster> ecalnecore_NonEM;
		std::vector<edm4hep::ConstCluster> ecalfrag;
		std::vector<edm4hep::ConstCluster> ecalundef;
		std::vector<edm4hep::ConstCluster> ecalfrag_TBM_CH;
		std::vector<edm4hep::ConstCluster> ecalfrag_TBM_NE_EM;
		std::vector<edm4hep::ConstCluster> ecalfrag_TBM_NE_NonEM;
		std::vector<edm4hep::ConstCluster> ecalundef_iso;
		std::vector<edm4hep::ConstCluster> ecalpotentialbackscattering;

		std::vector<edm4hep::ConstCluster> chargedclustercore;
		std::vector<edm4hep::ConstCluster> chargedclustercore_abs;

		std::vector<edm4hep::ConstCluster> selfmergedcluster; 
		std::vector<edm4hep::ConstCluster> non_chargedclustercore;
		std::vector<edm4hep::ConstCluster> onlyNeutralCore;

		std::vector<edm4hep::ConstCluster> non_charged_pem_neutral_core;
		std::vector<edm4hep::ConstCluster> pem_neutral_core;

		std::map<edm4hep::ConstTrack, int>MCPTrack_Type;
		std::map<edm4hep::ConstTrack, TVector3> Track_EndPoint;       //Last hit
		std::map<edm4hep::ConstTrack, TVector3> TrackStartPoint;
		std::map<edm4hep::ConstCluster, float> CluFD; 
		std::map<edm4hep::ConstCluster, float> CluT0;
		std::map<edm4hep::ConstCluster, float> Clu_Depth; 
		std::map<edm4hep::ConstCluster, TVector3> CluCoG;
	typedef DataHandle<edm4hep::MCParticleCollection> MCParticleColHandler;
	MCParticleColHandler m_mcParticle{"MCParticle", Gaudi::DataHandle::Reader, this};

	typedef DataHandle<edm4hep::TrackCollection> TrkColHandler;
	TrkColHandler m_trkcol {"MarlinTrkTracks", Gaudi::DataHandle::Reader, this};

	typedef DataHandle<edm4hep::ClusterCollection> CluColHandler;
	CluColHandler m_clucol {"EHBushes", Gaudi::DataHandle::Reader, this};

	DataHandle<edm4hep::CalorimeterHitCollection>   m_col_IsoHit{"IsoHits",Gaudi::DataHandle::Reader, this};

	DataHandle<edm4hep::ClusterCollection> m_1stclucol{"1stAbs",Gaudi::DataHandle::Writer, this};
	DataHandle<edm4hep::ClusterCollection> clucol_merged{"mergedCluC",Gaudi::DataHandle::Writer, this};
	DataHandle<edm4hep::ClusterCollection> m_mergedclu_chCol{"mergedCluCh",Gaudi::DataHandle::Writer, this};
	DataHandle<edm4hep::ClusterCollection> m_mergedclu_neCol{"mergedCluNe",Gaudi::DataHandle::Writer, this};
	DataHandle<edm4hep::ClusterCollection> m_chargedcoreclusterCol{"ChargedClu",Gaudi::DataHandle::Writer, this};
	DataHandle<edm4hep::ClusterCollection> NetralClu{"NetralClu",Gaudi::DataHandle::Writer, this};
	DataHandle<edm4hep::ReconstructedParticleCollection> m_chargeparticleCol{"ArborCharged",Gaudi::DataHandle::Writer, this};
	DataHandle<edm4hep::ReconstructedParticleCollection> nerecoparticleCol{"ArborNeutral",Gaudi::DataHandle::Writer, this};
	DataHandle<edm4hep::ReconstructedParticleCollection> m_arborrecoparticleCol{"ArborPFO",Gaudi::DataHandle::Writer, this};
     ArborToolLCIO * m_ArborToolLCIO;
	
};


#endif


