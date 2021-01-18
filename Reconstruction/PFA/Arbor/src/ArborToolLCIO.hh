#ifndef ARBORTOOLLCIO_H_
#define ARBORTOOLLCIO_H_


#include "k4FWCore/DataHandle.h"
#include "GaudiAlg/GaudiAlgorithm.h"
#include "edm4hep/ClusterCollection.h"
#include "edm4hep/ReconstructedParticleCollection.h"
#include "edm4hep/EventHeaderCollection.h"
#include "edm4hep/SimTrackerHitCollection.h"
#include "edm4hep/TrackerHitCollection.h"
#include "edm4hep/CalorimeterHitCollection.h"
#include "edm4hep/VertexCollection.h"
#include "edm4hep/TrackCollection.h"
#include "edm4hep/MCParticle.h" 
#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/MCRecoCaloAssociation.h"
#include "edm4hep/MCRecoTrackerAssociation.h"
#include "edm4hep/MCRecoTrackerAssociationCollection.h"
#include "edm4hep/MCRecoCaloAssociationCollection.h"
#include "edm4hep/MCRecoParticleAssociation.h"
#include "edm4hep/MCRecoParticleAssociationCollection.h"


#include <DDRec/DetectorData.h>
#include <DDRec/CellIDPositionConverter.h>
#include "DetInterface/IGeomSvc.h"

#include "TVector3.h"
#include "GaudiKernel/INTupleSvc.h"
#include "GaudiKernel/NTuple.h"
#include "GaudiKernel/MsgStream.h"
#include "ArborTool.h"
#include "HelixClassD.hh"

class CollectionMaps
{
public:
    CollectionMaps();
    void clear();
    std::map<std::string, std::vector<edm4hep::MCParticle> >     collectionMap_MC;
    std::map<std::string, std::vector<edm4hep::ConstCalorimeterHit> > collectionMap_CaloHit;
    std::map<std::string, std::vector<edm4hep::Vertex> >         collectionMap_Vertex;
    std::map<std::string, std::vector<edm4hep::Track> >          collectionMap_Track;
    std::map<std::string, std::vector<edm4hep::MCRecoCaloAssociation> > collectionMap_CaloRel;
    std::map<std::string, std::vector<edm4hep::MCRecoTrackerAssociation> > collectionMap_TrkRel;
};



class ArborToolLCIO  : public GaudiAlgorithm
{
	public:
		ArborToolLCIO(const std::string& name, ISvcLocator* svcLoc);
		//ArborToolLCIO(ISvcLocator* svcLoc);
		~ArborToolLCIO();
		typedef DataHandle<edm4hep::Cluster>  ClusterType;


		void Clean(); 
	void ClusterBuilding(DataHandle<edm4hep::ClusterCollection>& _currbranchcoll, std::vector<edm4hep::ConstCalorimeterHit> Hits, branchcoll BranchOrder, int DHCALFlag);
	
	float ClusterT0(edm4hep::ConstCluster a_Clu);

	int NHScaleV2( std::vector<edm4hep::ConstCalorimeterHit> clu0, int RatioX, int RatioY, int RatioZ );
	float FDV2( std::vector<edm4hep::ConstCalorimeterHit> clu);
	
	int NHScaleV3( edm4hep::ConstCluster clu0, int RatioX, int RatioY, int RatioZ );
	
	float FDV3( edm4hep::ConstCluster clu);
	
	int ActiveLayers(  std::vector<edm4hep::ConstCalorimeterHit> clu );
	
	float BushDis( edm4hep::ConstCluster clu1, edm4hep::ConstCluster clu2);
	
	
	float* SimpleDisTrackClu(edm4hep::ConstTrack a_trk, edm4hep::ConstCluster a_clu);
	
	float SimpleBushTimeTrackClu(edm4hep::ConstTrack a_trk, edm4hep::ConstCluster a_clu);
	
	int SimpleBushNC(edm4hep::ConstTrack a_trk, edm4hep::ConstCluster a_clu);
	
	float DisPointToBush( TVector3 Pos1, edm4hep::ConstCluster clu1);
	
	TVector3 ClusterCoG(edm4hep::ConstCluster inputCluser);
	
	edm4hep::ClusterCollection* ClusterVecMerge( std::vector<edm4hep::ConstCluster> inputClusters, TMatrixF ConnectorMatrix, DataHandle<edm4hep::ClusterCollection>& clucol );
	
	edm4hep::ClusterCollection* ClusterVecColl( std::vector<edm4hep::ConstCluster> inputClusters, DataHandle<edm4hep::ClusterCollection>& m_clucol);

	edm4hep::Cluster NaiveCluImpl(edm4hep::ConstCluster a0_clu);
	void NaiveCluConst(edm4hep::ConstCluster a0_clu, edm4hep::Cluster);
	
	std::vector<edm4hep::ConstCluster> CollClusterVec(const edm4hep::ClusterCollection * input_coll );
	
	std::vector<edm4hep::ConstCalorimeterHit> CollHitVec(const edm4hep::CalorimeterHitCollection * input_coll, float EnergyThreshold);
	
	std::vector<edm4hep::Cluster> ClusterHitAbsorbtion( std::vector<edm4hep::ConstCluster> MainClusters, std::vector<edm4hep::ConstCalorimeterHit> IsoHits, float DisThreshold );
	
	edm4hep::Cluster NaiveMergeClu(std::vector<edm4hep::ConstCluster> inputCluVec);

	void NaiveMergeCluConst(std::vector<edm4hep::ConstCluster> inputCluVec,edm4hep::Cluster MergedClu);
	std::vector<edm4hep::ConstCluster> ClusterAbsorbtion( std::vector<edm4hep::ConstCluster> MainClusters, std::vector<edm4hep::ConstCluster> FragClusters, float thresholds, float DepthSlope );
	
	
	int JointsBetweenBush(edm4hep::ConstCluster a_Clu, edm4hep::ConstCluster b_Clu, float CellSize);
	
	TVector3 CluEndP(edm4hep::ConstCluster a_clu);
	
	int ClusterFlag(edm4hep::ConstCluster a_tree, edm4hep::ConstTrack a_trk);
	
	int ClusterFlag1st(edm4hep::ConstCluster a_tree);
	
	int newPhotonTag(edm4hep::ConstCluster a_clu);
	float ClusterEE(edm4hep::ConstCluster inputCluster);
	
	float EMClusterEE( edm4hep::ConstCluster inputCluster );
	
	std::vector<float> ClusterTime(edm4hep::ConstCluster inputCluster);

	private:
	SmartIF<IGeomSvc> m_geosvc;//=service<IGeomSvc>("GeomSvc");
	dd4hep::DDSegmentation::BitFieldCoder*	m_decoder;// = m_geosvc->getDecoder("ECALBarrel");


//	m_geosvc= service<IGeomSvc>("GeomSvc");
//	m_decoder = m_geosvc->getDecoder("ECALBarrel");
	

};
#endif //
