#include "ArborToolLCIO.hh"
#include "ArborTool.h"
#include "Arbor.h"
#include "DetectorPos.hh"
#include "HelixClassD.hh"
#include <TMath.h>
#include <values.h>
#include <cmath>
#include <stdexcept>
#include <sstream>

#include "TVector3.h"
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

#include "DD4hep/Detector.h"
#include "DD4hep/IDDescriptor.h"
#include "DD4hep/Plugins.h"

#include <DDRec/DetectorData.h>
#include <DDRec/CellIDPositionConverter.h>
#include "DetInterface/IGeomSvc.h"

using namespace std;
/* 
void ClusterBuilding( LCEvent * evtPP, std::string Name, std::vector<CalorimeterHit*> Hits, std::vector< std::vector<int> > BranchOrder, int DHCALFlag )
{
	LCCollection *currbranchcoll = new LCCollectionVec(LCIO::CLUSTER);
	LCFlagImpl flag;
	flag.setBit(LCIO::CHBIT_LONG);
	currbranchcoll->setFlag(flag.getFlag());

	int NBranch = BranchOrder.size();
	int BranchSize = 0;
	float currBranchEnergy = 0;
	TVector3 SeedPos, currPos;
	float MinMag = 1E9;
	float currMag = 0; 
	float ECALTotalEn = 0; 
	float HCALTotalEn = 0;

	for(int i0 = 0; i0 < NBranch; i0++)
	{
		ClusterImpl* a_branch = new ClusterImpl();
		std::vector<int> currbranchorder = BranchOrder[i0];
		BranchSize = currbranchorder.size();
		currBranchEnergy = 0;
		ECALTotalEn = 0;
		HCALTotalEn = 0;
		// CalorimeterHit *Seedhit = Hits[currbranchorder[BranchSize - 1]];
		MinMag = 1E9;

		for(int j = 0; j < BranchSize; j++)
		{
			CalorimeterHit * a_hit = Hits[currbranchorder[j]];
			currPos = a_hit->getPosition();
			// currMag = DisSeedSurface(currPos);
			currMag = currPos.Mag();

			if( currMag < MinMag )
			{
				MinMag = currMag;
				SeedPos = currPos;
			}

			if(fabs(a_hit->getEnergy() - DHCALCalibrationConstant) < 1.0E-6 )
			{
				HCALTotalEn += DHCALCalibrationConstant;
			}
			else
			{
				ECALTotalEn += a_hit->getEnergy();
			}

			currBranchEnergy +=  a_hit->getEnergy();

			a_branch->addHit(a_hit, (float)1.0);
		}

		a_branch->setEnergy(currBranchEnergy);
		float ArraySeedPos[3] = { float(SeedPos.X()), float(SeedPos.Y()), float(SeedPos.Z()) };
		a_branch->setPosition( ArraySeedPos );

		a_branch->subdetectorEnergies().resize(6) ;
		a_branch->subdetectorEnergies()[0] = ECALTotalEn ;
		a_branch->subdetectorEnergies()[1] = HCALTotalEn ;
		a_branch->subdetectorEnergies()[2] = 0 ;
		a_branch->subdetectorEnergies()[3] = 0 ;
		a_branch->subdetectorEnergies()[4] = 0 ;
		a_branch->subdetectorEnergies()[5] = 0 ;

		currbranchcoll -> addElement(a_branch);
	}

	evtPP->addCollection(currbranchcoll, Name);
}
*/


ArborToolLCIO::ArborToolLCIO(const std::string& name,ISvcLocator* svcLoc)
     : GaudiAlgorithm(name, svcLoc)
{
	m_geosvc=service<IGeomSvc>("GeomSvc");
}
	

ArborToolLCIO::~ArborToolLCIO()
{
}
void ArborToolLCIO::ClusterBuilding( DataHandle<edm4hep::ClusterCollection>& _currbranchcoll, std::vector<edm4hep::ConstCalorimeterHit> Hits, branchcoll BranchOrder, int DHCALFlag )
{
	//DataHandle<edm4hep::ClusterCollection> _currbranchcoll {"Name",Gaudi::DataHandle::Writer, this};
	//DataHandle<edm4hep::ClusterCollection> _currbranchcoll=new ClusterType(Name, Gaudi::DataHandle::Writer, this);
	edm4hep::ClusterCollection* currbranchcoll = _currbranchcoll.createAndPut();

	int NBranch = BranchOrder.size();
	int BranchSize = 0;
	float currBranchEnergy = 0;
	TVector3 SeedPos, currPos;

	float MinMag = 1E9;
	float currMag = 0; 
	float ECALTotalEn = 0; 
	float HCALTotalEn = 0;

	for(int i0 = 0; i0 < NBranch; i0++)
	{
		auto a_branch=currbranchcoll->create();
		std::vector<int> currbranchorder = BranchOrder[i0];
		BranchSize = currbranchorder.size();
		currBranchEnergy = 0;
		ECALTotalEn = 0;
		HCALTotalEn = 0;
		MinMag = 1E9;

		for(int j = 0; j < BranchSize; j++)
		{
			auto a_hit= Hits[currbranchorder[j]];
			currPos =  TVector3(a_hit.getPosition().x,a_hit.getPosition().y,a_hit.getPosition().z);//a_hit.getPosition();
			currMag = currPos.Mag();

			if( currMag < MinMag )
			{
				MinMag = currMag;
				SeedPos = currPos;
			}

			if(fabs(a_hit.getEnergy() - DHCALCalibrationConstant) < 1.0E-6 )
			{
				HCALTotalEn += DHCALCalibrationConstant;
			}
			else
			{
				ECALTotalEn += a_hit.getEnergy();
			}

			currBranchEnergy +=  a_hit.getEnergy();

			a_branch.addToHits(a_hit);
		}

		a_branch.setEnergy(currBranchEnergy);
		float ArraySeedPos[3] = { float(SeedPos.X()), float(SeedPos.Y()), float(SeedPos.Z()) };
		a_branch.setPosition( ArraySeedPos );

		//a_branch.getSubdetectorEnergies().resize(6) ;
		//a_branch.getSubdetectorEnergies(0) = ECALTotalEn ;
		//a_branch.getSubdetectorEnergies(1) = HCALTotalEn ;
		//a_branch.getSubdetectorEnergies(2) = 0 ;
		//a_branch.getSubdetectorEnergies(3) = 0 ;
		//a_branch.getSubdetectorEnergies(4) = 0 ;
		//a_branch.getSubdetectorEnergies(5) = 0 ;

	}

}



int ArborToolLCIO::NHScaleV2( std::vector<edm4hep::ConstCalorimeterHit> clu0, int RatioX, int RatioY, int RatioZ )
{

	int ReScaledNH = 0;
	int NumHit = clu0.size();
	int tmpI = 0;
	int tmpJ = 0;
	int tmpK = 0;
	float tmpEn = 0;
	int NewCellID0 = 0;

	m_decoder = m_geosvc->getDecoder("EcalBarrelCollection");
	if(!m_decoder) m_decoder = m_geosvc->getDecoder("EcalEndcapsCollection");


	std::map <double, float> testIDtoEnergy;

	for(int i = 0; i < NumHit; i++)
	{
		auto hit = clu0[i];
		auto cellid = hit.getCellID();


		tmpI = m_decoder->get(cellid, "cellX")/RatioX;
		tmpJ = m_decoder->get(cellid, "cellY")/RatioY;
		tmpK = (m_decoder->get(cellid, "layer")+1)/RatioZ;
		tmpEn = hit.getEnergy();

		NewCellID0 = (tmpK<<24) + (tmpJ<<12) + tmpI;

		if(testIDtoEnergy.find(NewCellID0) == testIDtoEnergy.end() )
		{
			testIDtoEnergy[NewCellID0] = tmpEn;
		}
		else
		{
			testIDtoEnergy[NewCellID0] += tmpEn;
		}
	}

	ReScaledNH = testIDtoEnergy.size();

	return ReScaledNH;
}

float ArborToolLCIO::FDV2( std::vector<edm4hep::ConstCalorimeterHit> clu)
{
	float FractalDim = 0;
	int NReSizeHit[10] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
	int Scale[10] = {2, 3, 4, 5, 6, 7, 8, 9, 10, 20};
	int OriNHit = clu.size();

	for(int j = 0; j < 10; j++)
	{
		NReSizeHit[j] = NHScaleV2(clu, Scale[j], Scale[j], 1);
		FractalDim += 0.1 * TMath::Log(float(OriNHit)/NReSizeHit[j])/TMath::Log(float(Scale[j]));
	}

	if(clu.size() == 0) 
		FractalDim = -1; 

	return FractalDim;
}


int ArborToolLCIO::NHScaleV3( edm4hep::ConstCluster clu0, int RatioX, int RatioY, int RatioZ )
{

	int ReScaledNH = 0;
	int NumHit = clu0.hits_size();
	int tmpI = 0;
	int tmpJ = 0;
	int tmpK = 0;
	float tmpEn = 0;
	int NewCellID0 = 0;
	int NewCellID1 = 0;
	//m_geosvc=service<IGeomSvc>("GeomSvc");	

	m_decoder = m_geosvc->getDecoder("EcalBarrelCollection");
	if(!m_decoder) m_decoder = m_geosvc->getDecoder("EcalEndcapsCollection");

	std::map <double, float> testIDtoEnergy;
	double testlongID = 0;
	
	for(int i = 0; i < NumHit; i++)
	{
		auto hit = clu0.getHits(i);
		auto cellid = hit.getCellID();

		tmpI = m_decoder->get(cellid, "cellX")/RatioX;
		tmpJ = m_decoder->get(cellid, "cellY")/RatioY;
		tmpK = (m_decoder->get(cellid, "layer")+1)/RatioZ;
		tmpEn = hit.getEnergy();

		NewCellID0 = (tmpK<<24) + (tmpJ<<12) + tmpI;

		testlongID = NewCellID1*1073741824 + NewCellID0;
		if(testIDtoEnergy.find(testlongID) == testIDtoEnergy.end() )
		{
			testIDtoEnergy[testlongID] = tmpEn;
		}
		else
		{
			testIDtoEnergy[testlongID] += tmpEn;
		}
	}

	ReScaledNH = testIDtoEnergy.size();

	return ReScaledNH;

}

float ArborToolLCIO::FDV3( edm4hep::ConstCluster clu ){

	float FractalDim = -1;
        int NReSizeHit[10] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
	int Scale[5] = {2, 3, 4, 5, 6};
	int OriNHit = clu.hits_size();
	if(OriNHit > 0)
	{
		FractalDim = 0.0;
		for(int j = 0; j < 5; j++)
		{
			NReSizeHit[j] = NHScaleV3(clu, Scale[j], Scale[j], 1);
			FractalDim += 0.2 * TMath::Log(float(OriNHit)/NReSizeHit[j])/TMath::Log(float(Scale[j]));
		}
	}
	return FractalDim;
}


float ArborToolLCIO::BushDis( edm4hep::ConstCluster clu1, edm4hep::ConstCluster clu2)
{
	float DisBetweenBush = 1.0E10; 

	int cluSize1 = clu1.hits_size();
	int cluSize2 = clu2.hits_size();

	TVector3 HitPos1, HitPos2; 
	TVector3 PosDiff; 
	// TVector3 XXXPos; 

	for(int i = 0; i < cluSize1; i++)
	{
		HitPos1 = TVector3((clu1.getHits(i)).getPosition().x,(clu1.getHits(i)).getPosition().y,(clu1.getHits(i)).getPosition().z);
		for(int j = 0; j<cluSize2; j++)
		{
			HitPos2 = TVector3((clu2.getHits(j)).getPosition().x,(clu2.getHits(j)).getPosition().y,(clu2.getHits(j)).getPosition().z);
			PosDiff = HitPos1 - HitPos2;

			if(PosDiff.Mag() < DisBetweenBush )
			{
				DisBetweenBush = PosDiff.Mag();
			}
		}
	}


	return DisBetweenBush; 
}


float ArborToolLCIO::DisPointToBush(TVector3 Pos1, edm4hep::ConstCluster clu1)
{
	float Dis = 1.0E9; 
	float HitDis = 1.0E8;
	int clusize = clu1.hits_size();

	TVector3 HitPos; 

	for(int s = 0; s < clusize; s++)
	{
		HitPos = TVector3((clu1.getHits(s)).getPosition().x,(clu1.getHits(s)).getPosition().y,(clu1.getHits(s)).getPosition().z);
		HitDis = (HitPos - Pos1).Mag();
		if(HitDis < Dis) 
		{
			Dis = HitDis; 
		}
	}

	return Dis; 
}


TVector3 ArborToolLCIO::ClusterCoG(edm4hep::ConstCluster inputCluster)
{
	TVector3 CenterOfGravity; 

	int inputClusterSize = inputCluster.hits_size();

	TVector3 tmphitPos; 
	float tmphitEnergy;
	float sumhitEnergy = 0; 

	for(int i = 0; i < inputClusterSize; i++)
	{
		auto tmpHit = inputCluster.getHits(i);
		tmphitPos = TVector3(tmpHit.getPosition().x,tmpHit.getPosition().y,tmpHit.getPosition().z);
		tmphitEnergy = tmpHit.getEnergy();

		CenterOfGravity += tmphitPos*tmphitEnergy;
		sumhitEnergy += tmphitEnergy; 
	}

	CenterOfGravity = 1.0/sumhitEnergy * CenterOfGravity; 

	return CenterOfGravity; 
}


edm4hep::ClusterCollection* ArborToolLCIO::ClusterVecColl( std::vector<edm4hep::ConstCluster> inputClusters, DataHandle<edm4hep::ClusterCollection>& m_clucol )
{

	edm4hep::ClusterCollection* vec_coll_Clusters = m_clucol.createAndPut();

	int NClu = inputClusters.size();
	int CurrBranchSize = 0; 
	TVector3 SeedPos; 
	std::vector<float> CluEn; 
	std::vector<int> CluIndex; 

	for(int i0 = 0; i0 < NClu; i0++)
	{
		auto a_clu = inputClusters[i0];
		CluEn.push_back(a_clu.getEnergy());
	}
	CluIndex = SortMeasure(CluEn, 1);

	for(int i1 = 0; i1 < NClu; i1++)
	{
		auto branchtmp=vec_coll_Clusters->create();
		auto a_clu = inputClusters[CluIndex[i1]];

		CurrBranchSize = a_clu.hits_size();

		for(int j1 = 0; j1 < CurrBranchSize; j1 ++)
		{
			auto tmpHit = a_clu.getHits(j1);
			branchtmp.addToHits(tmpHit);
		}

		branchtmp.setPosition( a_clu.getPosition() );
		branchtmp.setEnergy(a_clu.getEnergy() );
		SeedPos = TVector3(a_clu.getPosition().x,a_clu.getPosition().y,a_clu.getPosition().z);
		branchtmp.setITheta( SeedPos.Theta() );       //To be replaced, those worse than 1st order appro
		branchtmp.setPhi( SeedPos.Phi()  );

	}

	return vec_coll_Clusters;
}

std::vector<edm4hep::ConstCluster> ArborToolLCIO::CollClusterVec(const edm4hep::ClusterCollection * input_coll )
{
	std::vector<edm4hep::ConstCluster> outputClusterVec; 


	outputClusterVec.clear();

	for(int i = 0; i < input_coll->size(); i++)	//We can have some sort here - according to depth/energy...
	{
		auto a_clu = (*input_coll)[i];
		outputClusterVec.push_back(a_clu);
	}

	return outputClusterVec; 
}


void ArborToolLCIO::NaiveCluConst(edm4hep::ConstCluster a0_clu,edm4hep::Cluster b0_clu)
{
	b0_clu.setPosition(a0_clu.getPosition());
	b0_clu.setEnergy(a0_clu.getEnergy());
	int NCaloHit = a0_clu.hits_size();
	float HitEn = 0; 
	float SubDEn[6] = {0, 0, 0, 0, 0, 0};

	for(int t0 = 0; t0 < NCaloHit; t0++)
	{
		auto a0_hit = a0_clu.getHits(t0);
		b0_clu.addToHits(a0_hit);
		HitEn = a0_hit.getEnergy();
		if(fabs(HitEn - DHCALCalibrationConstant) < 1.0E-6)
		{
			SubDEn[1] += HitEn;
		}
		else
		{	
			SubDEn[0] += HitEn;
		}
	}

	for(int i = 0; i < 6; i++)
	{
		b0_clu.addToSubdetectorEnergies(SubDEn[i]);
	}
	
}


edm4hep::Cluster ArborToolLCIO::NaiveCluImpl(edm4hep::ConstCluster a0_clu)
{
	edm4hep::Cluster b0_clu;
	b0_clu.setPosition(a0_clu.getPosition());
	b0_clu.setEnergy(a0_clu.getEnergy());
	int NCaloHit = a0_clu.hits_size();
	float HitEn = 0; 
	float SubDEn[6] = {0, 0, 0, 0, 0, 0};

	for(int t0 = 0; t0 < NCaloHit; t0++)
	{
		auto a0_hit = a0_clu.getHits(t0);
		b0_clu.addToHits(a0_hit);
		HitEn = a0_hit.getEnergy();
		if(fabs(HitEn - DHCALCalibrationConstant) < 1.0E-6)
		{
			SubDEn[1] += HitEn;
		}
		else
		{	
			SubDEn[0] += HitEn;
		}
	}

	for(int i = 0; i < 6; i++)
	{
		b0_clu.addToSubdetectorEnergies(SubDEn[i]);
	}
	
	return b0_clu; 
}

std::vector<edm4hep::ConstCalorimeterHit> ArborToolLCIO::CollHitVec(const edm4hep::CalorimeterHitCollection * input_coll, float EnergyThreshold)
{
	std::vector<edm4hep::ConstCalorimeterHit> outputHitVec;

	outputHitVec.clear();

	for(int i = 0; i < input_coll->size(); i++)      //We can have some sort here - according to depth/energy...
	{
		auto a_hit = (*input_coll)[i];
		if(a_hit.getEnergy() > EnergyThreshold)
		{
			outputHitVec.push_back(a_hit);
		}
	}

	return outputHitVec;
}


std::vector<edm4hep::Cluster> ArborToolLCIO::ClusterHitAbsorbtion( std::vector<edm4hep::ConstCluster> MainClusters, std::vector<edm4hep::ConstCalorimeterHit> IsoHits, float DisThreshold )	// Projective Distance + Hit Depth correlation; 
{
	std::vector<edm4hep::Cluster> outputClusterVec;

	int N_Core = MainClusters.size();
	int N_Hit = IsoHits.size();
	TVector3 HitPos, MBSeedPos;
	float currHitCoreDis = 0;  
	float MinHitCoreDis = 1.0E10; 
	int MinDisIndex = -1; 
	std::vector<std::pair<int, int> > Frag_Core_Links;
	std::pair<int, int> a_frag_core_link;

	for(int i0 = 0; i0 < N_Hit; i0++)
	{
		auto a_hit = IsoHits[i0];
		HitPos = TVector3(a_hit.getPosition().x,a_hit.getPosition().y,a_hit.getPosition().z);		
		MinHitCoreDis = 1.0E10;

		for(int j0 = 0; j0 < N_Core; j0++)
		{
			auto a_core = MainClusters[j0];
			currHitCoreDis = DisPointToBush(HitPos, a_core);
			if(currHitCoreDis < MinHitCoreDis)
			{
				MinHitCoreDis = currHitCoreDis; 
				MinDisIndex = j0;
			}
		}
		if(MinHitCoreDis < DisThreshold)
		{
			a_frag_core_link.first = i0;
			a_frag_core_link.second = MinDisIndex;
			Frag_Core_Links.push_back(a_frag_core_link);
		}
	}

	int N_frag_core_links = Frag_Core_Links.size();
	std::vector<edm4hep::ConstCalorimeterHit> tomerge_hits;
	float ClusterEn = 0; 

	for(int i2 = 0; i2 < N_Core; i2 ++)
	{
		auto a_core = MainClusters[i2];
		tomerge_hits.clear();

		for(int j4 = 0; j4 < N_frag_core_links; j4 ++)
		{
			a_frag_core_link = Frag_Core_Links[j4];
			if(a_frag_core_link.second == i2)
			{
				auto a_frag = IsoHits[a_frag_core_link.first];
				tomerge_hits.push_back(a_frag);
			}
		}
		edm4hep::Cluster a_mergedfrag_core;
		ClusterEn = 0; 

		for(unsigned int j2 = 0; j2 < a_core.hits_size(); j2++)
		{
			auto b_hit  = a_core.getHits(j2);
			a_mergedfrag_core.addToHits(b_hit);
			ClusterEn += b_hit.getEnergy();
		}

		for(unsigned int j3 = 0; j3 < tomerge_hits.size(); j3++)
		{
			auto c_hit = tomerge_hits[j3];
			a_mergedfrag_core.addToHits(c_hit);
			ClusterEn += c_hit.getEnergy();
		}

		a_mergedfrag_core.setPosition( a_core.getPosition() );
		a_mergedfrag_core.setEnergy(ClusterEn);
		MBSeedPos = TVector3(a_core.getPosition().x,a_core.getPosition().y,a_core.getPosition().z);
		a_mergedfrag_core.setITheta( MBSeedPos.Theta() );       //To be replaced, those worse than 1st order appro
		a_mergedfrag_core.setPhi( MBSeedPos.Phi()  );

		outputClusterVec.push_back(a_mergedfrag_core);
	}

	return outputClusterVec; 
}


void ArborToolLCIO::NaiveMergeCluConst(std::vector<edm4hep::ConstCluster> inputCluVec,edm4hep::Cluster MergedClu)
{

	int NClu = inputCluVec.size();
	float SeedDis = 1E9; 
	float MaxDis = 0; 
	float MergedCluEnergy = 0; 
	float HitEn = 0; 
	float SubDEn[6] = {0, 0, 0, 0, 0, 0};

	TVector3 CurrSeedPos, SeedPos, CurrHitPos, CluEndPos, CluRefDir;	//Seed Depth... CoG Comp...

	for(int i = 0; i < NClu; i++)
	{
		auto a_Clu = inputCluVec[i];
		CurrSeedPos =  TVector3(a_Clu.getPosition().x,a_Clu.getPosition().y,a_Clu.getPosition().z);
		MergedCluEnergy += a_Clu.getEnergy();

		if(CurrSeedPos.Mag() < SeedDis)
		{
			SeedPos = CurrSeedPos; 
			SeedDis = CurrSeedPos.Mag();
		}

		int CurrCluSize = a_Clu.hits_size();

		for(int j = 0; j < CurrCluSize; j++)
		{
			auto a_hit = a_Clu.getHits(j);
			MergedClu.addToHits(a_hit);
			CurrHitPos = TVector3(a_hit.getPosition().x,a_hit.getPosition().y,a_hit.getPosition().z);
			HitEn = a_hit.getEnergy();
			if(fabs(HitEn - DHCALCalibrationConstant) < 1.0E-6)	// ECAL, HCAL, Should use better criteria. 
			{
				SubDEn[1] += HitEn; 
			}
			else
			{
				SubDEn[0] += HitEn; 
			}

			if(CurrHitPos.Mag() > MaxDis)
			{
				MaxDis = CurrHitPos.Mag();
				CluEndPos = TVector3(a_hit.getPosition().x,a_hit.getPosition().y,a_hit.getPosition().z);	
			}
		}
	}
	CluRefDir = (CluEndPos - SeedPos);

	float ClusterSeedPos[3] = {float(SeedPos.X()), float(SeedPos.Y()), float(SeedPos.Z())};
	MergedClu.setPosition(ClusterSeedPos);
	MergedClu.setEnergy(MergedCluEnergy);
	MergedClu.setITheta( CluRefDir.Theta() );
	MergedClu.setPhi( CluRefDir.Phi() );

	//MergedClu.subdetectorEnergies().resize(6) ;
	for(int i = 0; i < 6; i++)
	{
	//	MergedClu.subdetectorEnergies()[i] = SubDEn[i];
		MergedClu.addToSubdetectorEnergies(SubDEn[i]);
	}

}
edm4hep::Cluster ArborToolLCIO::NaiveMergeClu(std::vector<edm4hep::ConstCluster> inputCluVec)
{
	edm4hep::Cluster MergedClu;

	int NClu = inputCluVec.size();
	float SeedDis = 1E9; 
	float MaxDis = 0; 
	float MergedCluEnergy = 0; 
	float HitEn = 0; 
	float SubDEn[6] = {0, 0, 0, 0, 0, 0};

	TVector3 CurrSeedPos, SeedPos, CurrHitPos, CluEndPos, CluRefDir;	//Seed Depth... CoG Comp...

	for(int i = 0; i < NClu; i++)
	{
		auto a_Clu = inputCluVec[i];
		CurrSeedPos =  TVector3(a_Clu.getPosition().x,a_Clu.getPosition().y,a_Clu.getPosition().z);
		MergedCluEnergy += a_Clu.getEnergy();

		if(CurrSeedPos.Mag() < SeedDis)
		{
			SeedPos = CurrSeedPos; 
			SeedDis = CurrSeedPos.Mag();
		}

		int CurrCluSize = a_Clu.hits_size();

		for(int j = 0; j < CurrCluSize; j++)
		{
			auto a_hit = a_Clu.getHits(j);
			MergedClu.addToHits(a_hit);
			CurrHitPos = TVector3(a_hit.getPosition().x,a_hit.getPosition().y,a_hit.getPosition().z);
			HitEn = a_hit.getEnergy();
			if(fabs(HitEn - DHCALCalibrationConstant) < 1.0E-6)	// ECAL, HCAL, Should use better criteria. 
			{
				SubDEn[1] += HitEn; 
			}
			else
			{
				SubDEn[0] += HitEn; 
			}

			if(CurrHitPos.Mag() > MaxDis)
			{
				MaxDis = CurrHitPos.Mag();
				CluEndPos = TVector3(a_hit.getPosition().x,a_hit.getPosition().y,a_hit.getPosition().z);	
			}
		}
	}
	CluRefDir = (CluEndPos - SeedPos);

	float ClusterSeedPos[3] = {float(SeedPos.X()), float(SeedPos.Y()), float(SeedPos.Z())};
	MergedClu.setPosition(ClusterSeedPos);
	MergedClu.setEnergy(MergedCluEnergy);
	MergedClu.setITheta( CluRefDir.Theta() );
	MergedClu.setPhi( CluRefDir.Phi() );

	//MergedClu.subdetectorEnergies().resize(6) ;
	for(int i = 0; i < 6; i++)
	{
	//	MergedClu.subdetectorEnergies()[i] = SubDEn[i];
		MergedClu.addToSubdetectorEnergies(SubDEn[i]);
	}

	return MergedClu;
}



std::vector<edm4hep::ConstCluster> ArborToolLCIO::ClusterAbsorbtion( std::vector<edm4hep::ConstCluster> MainClusters, std::vector<edm4hep::ConstCluster> FragClusters, float DisThreshold, float DepthSlope )	//ProjectiveDis
{
	std::vector<edm4hep::ConstCluster> outputClusterVec;

	int N_Core = MainClusters.size();
	int N_frag = FragClusters.size();

	//tag minimal distance

	float MinFragCoreDis = 1.0E10; 
	float CurrFragCoreDis = 0;
	int MinDisIndex = -1;  
	std::vector<std::pair<int, int> > Frag_Core_Links;
	std::pair<int, int> a_frag_core_link;
	std::map<int, int> TouchedFrag;
	TouchedFrag.clear();
	TVector3 fragPos; 

	for(int i0 = 0; i0 < N_frag; i0 ++)
	{
		auto a_frag = FragClusters[i0];
		fragPos = TVector3(a_frag.getPosition().x,a_frag.getPosition().y,a_frag.getPosition().z);
		MinFragCoreDis = 1.0E10;
		for(int j0 = 0; j0 < N_Core; j0++)
		{
			auto a_core = MainClusters[j0];
			CurrFragCoreDis = BushDis(a_frag, a_core);
			if(CurrFragCoreDis < MinFragCoreDis)
			{
				MinFragCoreDis = CurrFragCoreDis;
				MinDisIndex = j0;
			}
		}

		if( MinFragCoreDis < DisThreshold + DepthSlope*DisSeedSurface(fragPos))
		{
			a_frag_core_link.first = i0;
			a_frag_core_link.second = MinDisIndex;
			Frag_Core_Links.push_back(a_frag_core_link);
		}
	}

	int N_frag_core_links = Frag_Core_Links.size();
	std::vector<edm4hep::ConstCluster> tomerge_clu;

	for(int i4 = 0; i4 < N_Core; i4 ++)
	{
		auto a_core = MainClusters[i4];
		tomerge_clu.clear();
		tomerge_clu.push_back(a_core);

		for(int j4 = 0; j4 < N_frag_core_links; j4 ++)
		{
			a_frag_core_link = Frag_Core_Links[j4];
			if(a_frag_core_link.second == i4)
			{
				auto a_frag = FragClusters[a_frag_core_link.first];
				TouchedFrag[a_frag_core_link.first] = 1;
				tomerge_clu.push_back(a_frag);
			}
		}
		auto a_mergedfrag_core = NaiveMergeClu(tomerge_clu);
		edm4hep::ConstCluster a_mergedfrag_coreCon=a_mergedfrag_core;
		outputClusterVec.push_back(a_mergedfrag_core);
	}

	for(int i1 = 0; i1 < N_frag; i1++)
	{
		if(TouchedFrag.find(i1) == TouchedFrag.end())
		{
			auto a_frag = FragClusters[i1];
			//auto a_isofrag = NaiveCluImpl(a_frag);;
			outputClusterVec.push_back(a_frag);
		}
	}

	return outputClusterVec;
}


edm4hep::ClusterCollection* ArborToolLCIO::ClusterVecMerge( std::vector<edm4hep::ConstCluster> inputClusters, TMatrixF ConnectorMatrix, DataHandle<edm4hep::ClusterCollection>& clucol  )
{
	edm4hep::ClusterCollection* mergedbranches = clucol.createAndPut();


	int NinputClu = inputClusters.size();
	int Nrow = ConnectorMatrix.GetNrows();
	int Ncol = ConnectorMatrix.GetNcols();

	if(Ncol != NinputClu || Nrow != Ncol || Nrow != NinputClu)
	{
		cout<<"Size of Connector Matrix and inputClusterColl is not match"<<endl;
	}

	vector<edm4hep::ConstCluster> branchToMerge;
	edm4hep::ConstCluster Mergebranch_A;
	edm4hep::ConstCluster Mergebranch_B;
	edm4hep::ConstCluster tmpMergebranch;
	edm4hep::ConstCluster Mainbranch (0);

	TVector3 tmpClusterSeedPos, MBSeedPos;	

	int CurrBranchSize = 0;
	float SeedPosMin = 1.0E10;
	float BranchEnergy = 0;

	int FlagBranchTouch[Nrow];

	for(int i0 = 0; i0 < Nrow; i0++)
	{
		FlagBranchTouch[i0] = 0;
	}

	for(int ibran = 0; ibran < Nrow ; ibran++)
	{
		if(FlagBranchTouch[ibran] == 0)
		{
			Mergebranch_A = inputClusters[ibran];
			branchToMerge.push_back(Mergebranch_A);
			FlagBranchTouch[ibran] = 1;
			BranchEnergy = 0;
			auto branchtmp = mergedbranches->create();

			for(int jbran = ibran + 1; jbran < Nrow; jbran++)
			{
				if(FlagBranchTouch[jbran] == 0)
				{
					Mergebranch_B = inputClusters[jbran];
					if( ConnectorMatrix(ibran, jbran) > 0.1 )
					{
						branchToMerge.push_back(Mergebranch_B);
						FlagBranchTouch[jbran] = 1;
					}
				}
			}

			SeedPosMin = 1.0E10;

			for(unsigned int i1 = 0; i1 < branchToMerge.size(); i1++)
			{
				tmpMergebranch = branchToMerge[i1];
				tmpClusterSeedPos = TVector3(tmpMergebranch.getPosition().x,tmpMergebranch.getPosition().y,tmpMergebranch.getPosition().z);
				if( tmpClusterSeedPos.Mag() < SeedPosMin)
				{

					Mainbranch = tmpMergebranch;
					SeedPosMin = tmpClusterSeedPos.Mag();
				}

				CurrBranchSize = tmpMergebranch.hits_size();
				BranchEnergy += tmpMergebranch.getEnergy();

				for(int j1 = 0; j1 < CurrBranchSize; j1 ++)
				{
					auto tmpHit = tmpMergebranch.getHits(j1);
					branchtmp.addToHits(tmpHit);
				}
				
			}
			branchtmp.setPosition( Mainbranch.getPosition() );
			branchtmp.setEnergy(BranchEnergy);
			MBSeedPos = TVector3(Mainbranch.getPosition().x,Mainbranch.getPosition().y,Mainbranch.getPosition().z);
			branchtmp.setITheta( MBSeedPos.Theta() );       //To be replaced, those worse than 1st order appro
			branchtmp.setPhi( MBSeedPos.Phi()  );

			branchToMerge.clear();
		}
	}

	return mergedbranches;

}
int ArborToolLCIO::JointsBetweenBush(edm4hep::ConstCluster a_Clu, edm4hep::ConstCluster b_Clu, float CellSize)
{
	int NJoint = 0; 
	int a_CluSize = a_Clu.hits_size();
	int b_CluSize = b_Clu.hits_size();
	TVector3 aHitPos, bHitPos, PosDiff, aCluPos, bCluPos; 	
	aCluPos = TVector3(a_Clu.getPosition().x,a_Clu.getPosition().y,a_Clu.getPosition().z);
	bCluPos = TVector3(b_Clu.getPosition().x,b_Clu.getPosition().y,b_Clu.getPosition().z);

	for(int i = 0; i < a_CluSize; i++)
	{
		auto ahit = a_Clu.getHits(i);
		aHitPos = TVector3(ahit.getPosition().x,ahit.getPosition().y,ahit.getPosition().z);
		for(int j = 0; j < b_CluSize; j++)
		{
			auto bhit = b_Clu.getHits(j);
			bHitPos = TVector3(bhit.getPosition().x,bhit.getPosition().y,bhit.getPosition().z);
			PosDiff = aHitPos - bHitPos; 
			if(PosDiff.Mag() < 1.5*CellSize)	//allow Diag connect... else use 1.2
			{
				// if((aCluPos - bHitPos).Mag() < 60 || (bCluPos - aHitPos).Mag() < 60)
				NJoint++;	//Change to NJoint Hit...
			}
		}
	}

	return NJoint; 
}


float* ArborToolLCIO::SimpleDisTrackClu( edm4hep::ConstTrack a_trk, edm4hep::ConstCluster a_clu)
{
	float* Distance = new float[3];
       	Distance[0]	= 1.0E9;
       	Distance[1]	= 1.0E9;
       	Distance[2]	= 1.0E9;
	float minDis = 1.0E9;
	float BushDist[3] = {0, 0 ,0};
	HelixClassD * a_Helix = new HelixClassD();
	//float refPoint[3] = a_Helix->getReferencePoint();
	a_Helix->Initialize_Canonical(a_trk.getTrackStates(0).phi, a_trk.getTrackStates(0).D0, a_trk.getTrackStates(0).Z0, a_trk.getTrackStates(0).omega, a_trk.getTrackStates(0).tanLambda, BField);
	int NCaloHits = a_clu.hits_size();

	for(int i1 = 0; i1 < NCaloHits; i1++)
	{
		auto a_hit = a_clu.getHits(i1);

		float hitpos[3]={a_hit.getPosition().x,a_hit.getPosition().y,a_hit.getPosition().z};
		float BushTime = a_Helix->getDistanceToPoint(hitpos, BushDist);
		if(BushTime > 0 && BushDist[2] < minDis )
		{
			minDis = BushDist[2];
			Distance[0] = BushDist[0];
			Distance[1] = BushDist[1];
			Distance[2] = BushDist[2];
		}
	}
	delete a_Helix;

	return Distance;
}
float ArborToolLCIO::SimpleBushTimeTrackClu(edm4hep::ConstTrack a_trk, edm4hep::ConstCluster  a_clu)
{
        float Distance = 1.0E9;
        float Time = 0;
        float BushDist[3] = {0, 0 ,0};
        HelixClassD * a_Helix = new HelixClassD();
	a_Helix->Initialize_Canonical(a_trk.getTrackStates(0).phi, a_trk.getTrackStates(0).D0, a_trk.getTrackStates(0).Z0, a_trk.getTrackStates(0).omega, a_trk.getTrackStates(0).tanLambda, BField);
        int NCaloHits = a_clu.hits_size();

        for(int i1 = 0; i1 < NCaloHits; i1++)
        {
                auto a_hit = a_clu.getHits(i1);

		float hitpos[3]={a_hit.getPosition().x,a_hit.getPosition().y,a_hit.getPosition().z};
                float BushTime = a_Helix->getDistanceToPoint(hitpos, BushDist);
                if(BushTime > 0 && BushDist[2] < Distance )
                {
                        Time = BushTime;
			Distance = BushDist[2];
                }
        }
	delete a_Helix;
        return Time;
}

int ArborToolLCIO::SimpleBushNC(edm4hep::ConstTrack  a_trk, edm4hep::ConstCluster  a_clu)
{
	float Distance = 1.0E9;
	//float Time = 0; 
	int NC = 0;
	float BushDist[3] = {0, 0 ,0};
	HelixClassD * a_Helix = new HelixClassD();
	a_Helix->Initialize_Canonical(a_trk.getTrackStates(0).phi, a_trk.getTrackStates(0).D0, a_trk.getTrackStates(0).Z0, a_trk.getTrackStates(0).omega, a_trk.getTrackStates(0).tanLambda, BField);
	int NCaloHits = a_clu.hits_size();

	for(int i1 = 0; i1 < NCaloHits; i1++)
	{
		auto a_hit = a_clu.getHits(i1);

		float hitpos[3]={a_hit.getPosition().x,a_hit.getPosition().y,a_hit.getPosition().z};
		float BushTime = a_Helix->getDistanceToPoint(hitpos, BushDist);
		int nCircles = a_Helix->getNCircle(hitpos);
		if(BushTime > 0 && BushDist[2] < Distance )
		{
			//Time = BushTime;
			NC = nCircles;
			Distance = BushDist[2];
		}
	}

	delete a_Helix;
	return NC;
}

int ArborToolLCIO::ClusterFlag(edm4hep::ConstCluster a_tree, edm4hep::ConstTrack a_trk)
{
	// give each charged core cluster a flag
	//  Fragmentation:       999
	//  MIP: matched         130
	//       non-matched     131
	//  EM:  matched         110
	//       non-matched     111
	//  HAD: matched         211
	//       non-matched     212
	
	int CluIDFlag = 999;
	int ClusterID = 211;
	int EcalNHit, HcalNHit, NLEcal, NLHcal, NH[16];
	float avEnDisHtoL;
	float EcalEn, HcalEn, EClu, cluDepth, maxDepth, minDepth, MaxDisHel, MinDisHel, FD_all, FD_ECAL, FD_HCAL;
	float crdis, EEClu_L10, EEClu_R, EEClu_r, EEClu_p, rms_Ecal, rms_Hcal, rms_Ecal2, rms_Hcal2, av_NHE, av_NHH;
	int AL_Ecal, AL_Hcal;
	float FD_ECALF10;
	float dEdx = 0;

	const float mass = 0.139;       //Pion Mass

	HelixClassD * TrkInit_Helix = new HelixClassD();
	TrkInit_Helix->Initialize_Canonical(a_trk.getTrackStates(0).phi, a_trk.getTrackStates(0).D0, a_trk.getTrackStates(0).Z0, a_trk.getTrackStates(0).omega, a_trk.getTrackStates(0).tanLambda, BField);
	float TrackEn = mass*mass;


	for (int q3 = 0; q3 < 3; q3 ++)
	{
		TrackEn += (TrkInit_Helix->getMomentum()[q3])*(TrkInit_Helix->getMomentum()[q3]);
	}
	delete TrkInit_Helix;

	TrackEn = sqrt(TrackEn);
	int nSubTrk = a_trk.getTracks().size();

	//int NHit = a_tree->getHits().size();
	if(1 == 1) //if ( (NHit > 4 && TrackEn > 1) || TrackEn <= 1 )
	{

		TVector3 CluPos;
		CluPos = TVector3(a_tree.getPosition().x,a_tree.getPosition().y,a_tree.getPosition().z);
		TVector3 IntDir = ClusterCoG(a_tree)-CluPos;
		EClu = a_tree.getEnergy();
		EcalNHit = 0;
		HcalNHit = 0;
		EcalEn = 0;
		HcalEn = 0;
		float currDepth = 0;
		maxDepth = -100;
		minDepth = 1E6;
		MaxDisHel = -1;   //maximal distance from Track to Helix
		MinDisHel = 1E10;

		EEClu_R = 0;
		EEClu_r = 0;
		EEClu_p = 0;
		EEClu_L10 = 0;


		std::vector<edm4hep::ConstCalorimeterHit> Ecalhits;
		std::vector<edm4hep::ConstCalorimeterHit> Hcalhits;
		std::vector<edm4hep::ConstCalorimeterHit> allhits;
		std::vector<edm4hep::ConstCalorimeterHit> EH_1;
		std::vector<edm4hep::ConstCalorimeterHit> EH_2;
		std::vector<edm4hep::ConstCalorimeterHit> EH_3;
		std::vector<edm4hep::ConstCalorimeterHit> EH_4;
		std::vector<edm4hep::ConstCalorimeterHit> EH_5;
		std::vector<edm4hep::ConstCalorimeterHit> EH_6;
		std::vector<edm4hep::ConstCalorimeterHit> HH_1;
		std::vector<edm4hep::ConstCalorimeterHit> HH_2;
		std::vector<edm4hep::ConstCalorimeterHit> HH_3;
		std::vector<edm4hep::ConstCalorimeterHit> HH_4;
		std::vector<edm4hep::ConstCalorimeterHit> HH_5;
		std::vector<edm4hep::ConstCalorimeterHit> HH_6;
		std::vector<edm4hep::ConstCalorimeterHit> HH_7;
		std::vector<edm4hep::ConstCalorimeterHit> HH_8;
		std::vector<edm4hep::ConstCalorimeterHit> HH_9;
		std::vector<edm4hep::ConstCalorimeterHit> HH_0;
		std::vector<edm4hep::ConstCalorimeterHit> Ecalf10hits;



		allhits.clear();
		Ecalhits.clear();
		Hcalhits.clear();
		EH_1.clear();
		EH_2.clear();
		EH_3.clear();
		EH_4.clear();
		EH_5.clear();
		EH_6.clear();
		HH_1.clear();
		HH_2.clear();
		HH_3.clear();
		HH_4.clear();
		HH_5.clear();
		HH_6.clear();
		HH_7.clear();
		HH_8.clear();
		HH_9.clear();
		HH_0.clear();
		Ecalf10hits.clear();


		HelixClassD * currHelix = new HelixClassD();
		currHelix->Initialize_Canonical(a_trk.getTrackStates(0).phi, a_trk.getTrackStates(0).D0, a_trk.getTrackStates(0).Z0, a_trk.getTrackStates(0).omega, a_trk.getTrackStates(0).tanLambda, BField);
		float BushDist[3] = {0, 0, 0};
		float BushTime = 0;

		std::vector<float> hitTheta;
		hitTheta.clear();

		for(unsigned int j1 = 0; j1 < a_tree.hits_size(); j1++)
		{
			auto a_hit = a_tree.getHits(j1);
			float hitpos[3]={a_hit.getPosition().x,a_hit.getPosition().y,a_hit.getPosition().z};
			BushTime = currHelix->getDistanceToPoint(hitpos, BushDist);
			TVector3 tmpPos = TVector3(a_hit.getPosition().x,a_hit.getPosition().y,a_hit.getPosition().z);
			hitTheta.push_back(tmpPos.Theta());
			if(BushTime > 0)
			{
				if(BushDist[2] > MaxDisHel)
				{
					MaxDisHel = BushDist[2];
				}
				if(BushDist[2] < MinDisHel)
				{
					MinDisHel = BushDist[2];
				}
			}
		}
		delete currHelix;

		float totTheta = 0;
		float avTheta = 0;
		float SDTheta;

		for(int t0 = 0; t0 < int(hitTheta.size()); t0++)
		{
			float tmpTheta = hitTheta[t0];
			totTheta += tmpTheta;
		}

		avTheta = totTheta/float(hitTheta.size());
		SDTheta = 0;

		for(int t1 = 0; t1 < int(hitTheta.size()); t1++)
		{
			float tmpTheta = hitTheta[t1];
			SDTheta += pow(tmpTheta-avTheta,2);
		}
		SDTheta = sqrt(SDTheta/float(hitTheta.size()));

		TVector3 HitPos;
		int currCluNHits = a_tree.hits_size();
		if(currCluNHits == 0) return CluIDFlag;
		int index1 = 0, index2 = 0;

		for(int s1 = 0; s1 < currCluNHits; s1++)
		{
			auto a_hit = a_tree.getHits(s1);
			allhits.push_back(a_hit);
			auto cellid= a_hit.getCellID();
			int NLayer =  m_decoder->get(cellid, "layer");

			HitPos = TVector3(a_hit.getPosition().x,a_hit.getPosition().y,a_hit.getPosition().z);

			currDepth = DisSeedSurface(HitPos);
			crdis = (CluPos-HitPos).Mag()*sin((CluPos-HitPos).Angle(IntDir));


			if(currDepth > maxDepth)
			{
				maxDepth = currDepth;
				index1 = s1;
			}
			if(currDepth < minDepth)
			{
				minDepth = currDepth;
				index2 = s1;
			}

			if( fabs(a_hit.getEnergy() - DHCALCalibrationConstant) < 1.0E-6 )      //or other fancy judgements...^M
			{
				HcalNHit++;
				HcalEn += a_hit.getEnergy();
				Hcalhits.push_back(a_hit);
				if(NLayer < 5)
				{
					HH_1.push_back(a_hit);
				}
				else if(NLayer < 10)
				{
					HH_2.push_back(a_hit);
				}
				else if(NLayer < 15)
				{
					HH_3.push_back(a_hit);
				}
				else if(NLayer < 20)
				{
					HH_4.push_back(a_hit);
				}
				else if(NLayer < 25)
				{
					HH_5.push_back(a_hit);
				}
				else if(NLayer < 30)
				{
					HH_6.push_back(a_hit);
				}
				else if(NLayer < 35)
				{
					HH_7.push_back(a_hit);
				}
				else if(NLayer < 40)
				{
					HH_8.push_back(a_hit);
				}
				else if(NLayer < 45)
				{
					HH_9.push_back(a_hit);
				}
				else
				{
					HH_0.push_back(a_hit);
				}
			}
			else
			{
				EcalNHit++;
				EcalEn += a_hit.getEnergy();
				Ecalhits.push_back(a_hit);
				if(NLayer< 10) Ecalf10hits.push_back(a_hit);
				if(crdis < 22) EEClu_R += a_hit.getEnergy();
				if(crdis < 11) EEClu_r += a_hit.getEnergy();
				if(crdis < 6) EEClu_p += a_hit.getEnergy();
				if(NLayer < 5)
				{
					EH_1.push_back(a_hit);
					EEClu_L10 += a_hit.getEnergy();
				}
				else if(NLayer < 10)
				{
					EH_2.push_back(a_hit);
					EEClu_L10 += a_hit.getEnergy();
				}
				else if(NLayer < 15)
				{
					EH_3.push_back(a_hit);
				} 
				else if(NLayer < 20)
				{
					EH_4.push_back(a_hit);
				} 
				else if(NLayer < 25)
				{
					EH_5.push_back(a_hit);
				}
				else
				{
					EH_6.push_back(a_hit);
				}
			}
		}

		auto maxdis_hit = a_tree.getHits(index1);
		auto mindis_hit = a_tree.getHits(index2);
		TVector3 maxpos = TVector3(maxdis_hit.getPosition().x,maxdis_hit.getPosition().y,maxdis_hit.getPosition().z);
		TVector3 minpos = TVector3(mindis_hit.getPosition().x,mindis_hit.getPosition().y,mindis_hit.getPosition().z);
		TVector3 GraPos = ClusterCoG(a_tree);
		cluDepth = (maxpos-minpos).Mag();

		float totHitEn = 0;
		float totHitEnDis = 0;
		float HitEn;

		for(int s2 = 0; s2 < currCluNHits; s2++)
		{
			auto a_hit2 = a_tree.getHits(s2);
			HitPos = TVector3(a_hit2.getPosition().x,a_hit2.getPosition().y,a_hit2.getPosition().z);
			HitEn  = a_hit2.getEnergy();
			TVector3 par1 = GraPos-minpos;
			TVector3 par2 = minpos-HitPos;
			TVector3 par3 = par1.Cross(par2);
			float disHtoL = par3.Mag()/par1.Mag();
			totHitEn+=HitEn;
			totHitEnDis+=HitEn*disHtoL;
		}
		avEnDisHtoL = totHitEnDis/totHitEn;
		FD_all = FDV2(allhits);
		FD_ECAL = FDV2(Ecalhits);
		FD_HCAL = FDV2(Hcalhits);
		FD_ECALF10 = FDV2(Ecalf10hits);

		NLEcal = 0;
		NLHcal = 0;
		for(int p0 = 0; p0 < 8; p0++)
		{
			NH[p0] = 0;
		}

		NH[0] = EH_1.size();
		NH[1] = EH_2.size();
		NH[2] = EH_3.size();
		NH[3] = EH_4.size();
		NH[4] = EH_5.size();
		NH[5] = EH_6.size();
		NH[6] = HH_1.size();
		NH[7] = HH_2.size();
		NH[8] = HH_3.size();
		NH[9] = HH_4.size();
		NH[10] = HH_5.size();
		NH[11] = HH_6.size();
		NH[12] = HH_7.size();
		NH[13] = HH_8.size();
		NH[14] = HH_9.size();
		NH[15] = HH_0.size();

		NLEcal = ActiveLayers(Ecalhits);
		NLHcal = ActiveLayers(Hcalhits);
		cout<<"NLEcal "<<NLEcal<<" "<<Ecalhits.size()<<" "<<endl;


		float sum_NHE = 0, sum_NHH = 0;
		av_NHE = 0;
		av_NHH = 0;
		AL_Ecal = 0;
		AL_Hcal = 0;

		for(int r1 = 0; r1 < 16; r1++)
		{
			if(r1 < 6 && NH[r1]>0)
			{
				sum_NHE += NH[r1];
				AL_Ecal++;
			}
			if(r1 >= 6 && NH[r1]>0)
			{
				sum_NHH += NH[r1];
				AL_Hcal++;
			}
		}
		if(AL_Ecal > 0)
			av_NHE = sum_NHE/AL_Ecal;
		if(AL_Hcal > 0)
			av_NHH = sum_NHH/AL_Hcal;

		rms_Ecal = 0;
		rms_Hcal = 0;
		rms_Ecal2 = 0;
		rms_Hcal2 = 0;	
		for(int r0 = 0; r0 < 16; r0++)
		{
			if(r0 < 6)
			{
				if(NH[r0] > 0)
				{
					rms_Ecal+=pow(NH[r0]-av_NHE,2);
					rms_Ecal2 += pow(NH[r0],2);
				}
			}
			else
			{
				if(NH[r0] > 0)
				{
					rms_Hcal+=pow(NH[r0]-av_NHH,2);
					rms_Hcal2 += pow(NH[r0],2);
				}
			}
		}
		if(AL_Ecal > 0)
		{
			rms_Ecal2 = sqrt(rms_Ecal2/AL_Ecal);
			rms_Ecal = sqrt(rms_Ecal/AL_Ecal);
		}
		else
		{
			rms_Ecal2 = -1;
			rms_Ecal = -1;
		}
		if(AL_Hcal > 0)
		{
			rms_Hcal2 = sqrt(rms_Hcal2/AL_Ecal);
			rms_Hcal = sqrt(rms_Hcal/AL_Ecal);
		}
		else
		{
			rms_Hcal2 = -1;
			rms_Hcal = -1;
		}




		bool cutmu1 = EcalNHit+2.*HcalNHit < 500;
		bool cutmu2;

		if(cluDepth < 1100) cutmu2 = 0;
		else if(cluDepth < 1400) cutmu2 = FD_all < 0.5/pow(400,2)*pow((cluDepth-1000),2);
		else cutmu2 = 1;

		bool cutmu3;
		if(TrackEn > 70) cutmu3 = FD_all < (600./avEnDisHtoL + 20)/100.;
		else if (TrackEn > 10) cutmu3 = FD_all < (600./avEnDisHtoL -10 + 0.5*TrackEn)/100.;
		else if (TrackEn > 7.5) cutmu3 = FD_all < (600./avEnDisHtoL -10)/100.;
		else cutmu3 = 1;
		bool cutmu3b;
		if (TrackEn > 10) cutmu3b = avEnDisHtoL < 25;
		else if (TrackEn > 4.5) cutmu3b = avEnDisHtoL < 25 + 10 * (10 - TrackEn);
		else cutmu3b = 1;
		bool cutmu10en = cutmu1 && cutmu2 && cutmu3 && cutmu3b;


		bool cutmu4 = FD_HCAL >= 0;
		bool cutmu5 = cluDepth > 750 - 20/TrackEn;
		bool cutmu6 = cluDepth > 1200;
		bool cutmu7 = FD_all < 0.3/sqrt(400.)*sqrt(cluDepth-1200.);
		bool cutmu8 = rms_Hcal < 10 && FD_HCAL < -0.25/600.*MaxDisHel+0.25;
		bool cutmu9 = FD_all < 0.35/sqrt(400)*sqrt(400-MaxDisHel) && MaxDisHel < 400 && EcalNHit+2.*HcalNHit > 85;
		bool cutmu;

		if(TrackEn > 9.5) cutmu = cutmu10en;
		else if(TrackEn < 1.5)
		{
			if(cutmu5) cutmu = 1;
			else cutmu = 0;
		}
		else if(TrackEn < 3.5)
		{
			if(cutmu4 && cutmu6 && cutmu7) cutmu = 1;
			else cutmu = 0;
		}
		else
		{
			if(cutmu3b && cutmu4 && cutmu6 && cutmu7 && cutmu8 && cutmu9) cutmu = 1;
			else cutmu = 0;
		}

		//bool cute1 = EcalEn/(EcalEn+HcalEn) > 0.9;
		bool cute2 = (dEdx > 0.17e-6 && dEdx < 0.3e-6)||dEdx==0;
		bool cute3 = FD_all > 0.9*sin(cluDepth*1.57/800.) && cluDepth < 800;
		bool cute4 = FD_ECALF10 > 0.9*sin((cluDepth-200)*1.57/600.);
		bool cute5 = EEClu_r/TrackEn > 0.8 * sin((avEnDisHtoL-15)*1.57/20.) && avEnDisHtoL < 35;
		bool cute6 = FD_ECAL > 0.2 * log10(EcalNHit) && log10(EcalNHit) > 1.5;
		bool cute7 = FD_all >= 0.6/1.5*(log10(EcalNHit+2.*HcalNHit)-1.2-0.4*TrackEn/100.);
		bool cute8 = SDTheta < 0.012/1.5*(log10(EcalNHit+2.*HcalNHit)-1);
		bool cute;
		if(TrackEn < 1.5) cute = cute2;
		else cute = cute2 && cute3 && cute4 && cute5 && cute6 && cute7 && cute8;


		if(cutmu) ClusterID = 13;
		else if (cute)  ClusterID = 11;

		if(ClusterID == 13 )
		{
			if(NLEcal+NLHcal < 30) CluIDFlag = 131;
			else CluIDFlag = 130;
		}

		if(ClusterID == 11 )
		{
			if(EClu < 0.8*TrackEn) CluIDFlag = 111;
			else CluIDFlag = 110;
		}

		if(ClusterID == 211 )
		{
			if(EClu < 0.8*TrackEn) CluIDFlag = 212;
			else CluIDFlag = 211;
		}

	}

	return CluIDFlag;

}


int ArborToolLCIO::ActiveLayers(  std::vector<edm4hep::ConstCalorimeterHit> clu )
{
	std::vector<int> hitlayers; 
	hitlayers.clear();

	int NHits = clu.size();	
	int tmpK = 0;	//Layer Number
	int tmpS = 0; 
	int tmpID = 0;
	
	for(int i = 0; i < NHits; i++)
	{
		auto hit = clu[i];
		auto cellid= hit.getCellID();
		tmpK = m_decoder->get(cellid, "layer")+1 ;
		tmpS = m_decoder->get(cellid, "stave")+1 ;
		// cout<<"tmpK "<<tmpK<<endl; 
		tmpID = tmpS * 50 + tmpK;

		if( std::find(hitlayers.begin(), hitlayers.end(), tmpID) == hitlayers.end() )
		{
			hitlayers.push_back(tmpID);
		}
	}

	return hitlayers.size();
}


float ArborToolLCIO::ClusterT0(edm4hep::ConstCluster a_Clu)
{
	float T0 = 1E9; 
	float tmpTime = 0; 
	TVector3 CluHitPos; 
	for(unsigned int i = 0; i < a_Clu.hits_size(); i++)
	{
		auto a_hit = a_Clu.getHits(i);
		CluHitPos = TVector3(a_hit.getPosition().x,a_hit.getPosition().y,a_hit.getPosition().z);

		tmpTime = a_hit.getTime() - 1.0/300*CluHitPos.Mag();
		if(tmpTime < T0 && tmpTime > 0)
		{
			T0 = tmpTime; 
		}
	}
	return T0;
}


int ArborToolLCIO::newPhotonTag(edm4hep::ConstCluster a_clu)
{
	int flag=0;

	TVector3 PosClu = TVector3(a_clu.getPosition().x,a_clu.getPosition().y,a_clu.getPosition().z);
	float Depth = 0; 
	Depth = DisSeedSurface(PosClu);

	int CluSize= a_clu.hits_size();
	float ClusFD=FDV3(a_clu);
	float ClusT0=ClusterT0(a_clu);

	
	if(ClusFD>0.18*TMath::Log((float)CluSize)-0.53&&ClusFD<0.16*TMath::Log((float)CluSize)+0.025&&ClusFD>-0.2*TMath::Log((float)CluSize)+0.4&&((log10(ClusT0)<-2&&log10(Depth)<2&&log10(CluSize)>2)||(log10(ClusT0)<-1.5&&log10(CluSize)<2)))
	{
		flag=1;
	}
	

	return flag;



}


int ArborToolLCIO::ClusterFlag1st(edm4hep::ConstCluster a_tree)
{
	int ClusterID = 211;
	int EcalNHit, HcalNHit, NH_ECALF10;
	float avEnDisHtoL;
	float EcalEn, HcalEn, EClu, cluDepth, maxDepth, minDepth, FD_all, FD_ECAL, FD_HCAL;
	float FD_ECALF10;
	

	TVector3 CluPos;
	CluPos = TVector3(a_tree.getPosition().x,a_tree.getPosition().y,a_tree.getPosition().z);

	EClu = a_tree.getEnergy();
	EcalNHit = 0;
	HcalNHit = 0;
	EcalEn = 0;
	HcalEn = 0;
	float currDepth = 0;
	maxDepth = -100;
	minDepth = 1E6;

	std::vector<edm4hep::ConstCalorimeterHit> allhits;
	std::vector<edm4hep::ConstCalorimeterHit> Ecalhits;
	std::vector<edm4hep::ConstCalorimeterHit> Hcalhits;
	std::vector<edm4hep::ConstCalorimeterHit> Ecalf10hits;

	allhits.clear();
	Ecalhits.clear();
	Hcalhits.clear();
	Ecalf10hits.clear();

	std::vector<float> hitTheta;
	hitTheta.clear();

	for(unsigned int j1 = 0; j1 < a_tree.hits_size(); j1++)
	{
		auto a_hit = a_tree.getHits(j1);

		TVector3 tmpPos = TVector3(a_hit.getPosition().x,a_hit.getPosition().y,a_hit.getPosition().z);
		hitTheta.push_back(tmpPos.Theta());

	}

	float totTheta = 0;
	float avTheta = 0;
	float SDTheta;

	for(int t0 = 0; t0 < int(hitTheta.size()); t0++)
	{
		float tmpTheta = hitTheta[t0];
		totTheta += tmpTheta;
	}

	avTheta = totTheta/float(hitTheta.size());
	SDTheta = 0;

	for(int t1 = 0; t1 < int(hitTheta.size()); t1++)
	{
		float tmpTheta = hitTheta[t1];
		SDTheta += pow(tmpTheta-avTheta,2);
	}
	SDTheta = sqrt(SDTheta/float(hitTheta.size()));

	TVector3 HitPos;
	int currCluNHits = a_tree.hits_size();
	if(currCluNHits == 0) return 1;
	int index1 = 0, index2 = 0;
	int HitDepth50 = 1.0e4;

	for(int s1 = 0; s1 < currCluNHits; s1++)
	{
		auto a_hit = a_tree.getHits(s1);
		allhits.push_back(a_hit);
		auto cellid = a_hit.getCellID();
		int NLayer = m_decoder->get(cellid, "layer");
		float tmpHitEn = a_hit.getEnergy();
		HitPos = TVector3(a_hit.getPosition().x,a_hit.getPosition().y,a_hit.getPosition().z);

		currDepth = DisSeedSurface(HitPos);
		if(tmpHitEn > 0.05 && currDepth < HitDepth50)
		{
			HitDepth50 = currDepth;
		}


		if(currDepth > maxDepth)
		{
			maxDepth = currDepth;
			index1 = s1;
		}
		if(currDepth < minDepth)
		{
			minDepth = currDepth;
			index2 = s1;
		}

		if( fabs(a_hit.getEnergy() - DHCALCalibrationConstant) < 1.0E-6 )      //or other fancy judgements...^M
		{
			HcalNHit++;
			HcalEn += a_hit.getEnergy();
			Hcalhits.push_back(a_hit);
		}
		else
		{
			EcalNHit++;
			EcalEn += a_hit.getEnergy();
			Ecalhits.push_back(a_hit);
			if(NLayer< 10) Ecalf10hits.push_back(a_hit);
		}
	}
	auto maxdis_hit = a_tree.getHits(index1);
	auto mindis_hit = a_tree.getHits(index2);

		TVector3 maxpos = TVector3(maxdis_hit.getPosition().x,maxdis_hit.getPosition().y,maxdis_hit.getPosition().z);
		TVector3 minpos = TVector3(mindis_hit.getPosition().x,mindis_hit.getPosition().y,mindis_hit.getPosition().z);
	TVector3 GraPos = ClusterCoG(a_tree);
	cluDepth = (maxpos-minpos).Mag();

	float totHitEn = 0;
	float totHitEnDis = 0;
	float HitEn;

	for(int s2 = 0; s2 < currCluNHits; s2++)
	{
		auto  a_hit2 = a_tree.getHits(s2);
		HitPos = TVector3(a_hit2.getPosition().x,a_hit2.getPosition().y,a_hit2.getPosition().z);
		HitEn  = a_hit2.getEnergy();
		TVector3 par1 = GraPos-minpos;
		TVector3 par2 = minpos-HitPos;
		TVector3 par3 = par1.Cross(par2);
		float disHtoL = par3.Mag()/par1.Mag();
		totHitEn+=HitEn;
		totHitEnDis+=HitEn*disHtoL;
	}
	avEnDisHtoL = totHitEnDis/totHitEn;
	FD_all = FDV2(allhits);
	FD_ECAL = FDV2(Ecalhits);
	FD_HCAL = FDV2(Hcalhits);
	FD_ECALF10 = FDV2(Ecalf10hits);
	NH_ECALF10 = Ecalf10hits.size();


	bool cute1 = log10(FD_all/log10(EClu)) > (log10(cluDepth/log10(EClu))-2.9);
	bool cute2 = FD_ECAL > 0.1 && FD_ECAL > 0.2*log10(avEnDisHtoL-3)+0.05*sqrt(avEnDisHtoL-3);
	float x;
	if (EcalNHit == 0) x = 0;
	else x = float(NH_ECALF10)/float(EcalNHit);
	bool cute3;
	if (x < 0.9)  cute3 = FD_ECALF10 > 0.3/sqrt(sqrt(0.9))*sqrt(sqrt(0.9-x));
	else cute3 = 1;
	bool cute4 = HcalNHit/EClu < 0.3;                
	bool cute;
	cute = cute1 && cute2 && cute3 && cute4;

	bool cutmu1 = (cluDepth > 1300 || EClu < 4/1000.*cluDepth+0.9);
	bool cutmu2 = EClu < 15;
	bool cutmu3 = avEnDisHtoL*SDTheta/EClu <0.02 && FD_all < 0.4/0.025*(0.025-avEnDisHtoL*SDTheta/EClu);
	bool cutmu;

	cutmu = cutmu1 && cutmu2 && cutmu3;


	if(currCluNHits <= 4) ClusterID = 1;
	else if(cute) ClusterID = 11;
	else if (cutmu)  ClusterID = 13;
	else 
	{
		bool cutef1 = FD_HCAL == -1;
		bool cutef2 = minDepth < 50;
		bool cutef3 = FD_ECAL > 0.2;
		bool cutef3b = 1;
		if (FD_ECAL < 0.6) cutef3b = FD_ECAL > 6/3.5*(1.75-log10((EcalNHit+HcalNHit)/EClu));
		bool cutef4;
		if (avEnDisHtoL <= 8)
			cutef4 = FD_ECALF10 > 0.2;
		else
			cutef4 = FD_ECALF10 > 0.2+0.8*sqrt(log10(avEnDisHtoL/8));
		bool cutef5 = FD_ECAL/EClu > (log10(EcalNHit)-1.6)*(log10(EcalNHit)-1.6)*4+0.25 && FD_ECALF10 > (avEnDisHtoL-9)*(avEnDisHtoL-9)*0.006+0.2;
		bool cutef;

		cutef = (cutef1 && cutef2 && cutef3 && cutef3b && cutef4) || cutef5;
		if(cutef) ClusterID = 11;  
	}

	if(ClusterID == 211)
	{
		bool cutmuf1 = FD_all == 0 && log10((EcalNHit+2*HcalNHit)*EClu) > 1.2;
		bool cutmuf2 = log10((EcalNHit+2*HcalNHit)*EClu) > 1.9 && log10((EcalNHit+2*HcalNHit)*EClu) < 3.1;
		bool cutmuf3 = log10(FD_all/cluDepth) < -223.074+431.315*log10((EcalNHit+2*HcalNHit)*EClu)-331.72*pow(log10((EcalNHit+2*HcalNHit)*EClu),2)+124.588*pow(log10((EcalNHit+2*HcalNHit)*EClu),3)-22.8786*pow(log10((EcalNHit+2*HcalNHit)*EClu),4)+1.64577*pow(log10((EcalNHit+2*HcalNHit)*EClu),5);
		bool cutmuf = cutmuf1 || (cutmuf2 && cutmuf3);
		if(cutmuf) ClusterID = 13;
		else if(minDepth < 0.77+0.23*EClu) ClusterID = 211;
		else ClusterID = 2; 
	}

	return ClusterID;
}

float ArborToolLCIO::ClusterEE(edm4hep::ConstCluster inputCluster)
{
	float ClusterEnergy = 0;
	float tmpCluEn = 0;
	int CluType = 211;
	if(ClusterFlag1st(inputCluster) == 11) CluType = 22;

	int NCluHit = inputCluster.hits_size();
	float hitEn = 0;
	float EnCorrector = 1;

	if(CluType == 22)
	{
		ClusterEnergy = EMClusterEE(inputCluster);
	}
	else
	{
		for(int i = 0; i < NCluHit; i++)
		{
			auto a_hit = inputCluster.getHits(i);
			hitEn = a_hit.getEnergy();

			if(CluType == 211)
			{
				if( hitEn < 1.5 )       // Or some function of the initCluEn
				{
					tmpCluEn += hitEn;
				}
				else
				{
					cout<<"Ultra Hot Hit"<<endl;	// Use ShuZhen's Hot Hit Finding function...
					tmpCluEn += 0.05;       // MIP Hit Energy, Value to be determined
				}
			}
			else    // For EM & MIP: should veto accordingly
			{
				tmpCluEn += hitEn;
			}
		}
		
		EnCorrector = 1;
		if(tmpCluEn > 1.5 && tmpCluEn < 22 && 1 == 0)
		{
			EnCorrector = 0.6*(1 + 1.0/log10(tmpCluEn)); 
		}
		ClusterEnergy = tmpCluEn*EnCorrector;
		//ClusterEnergy = HADClusterEE(tmpCluEn,inputCluster);
	}

	return ClusterEnergy;
}


float ArborToolLCIO::EMClusterEE( edm4hep::ConstCluster inputCluster )
{
	float aaa = -50.2288, ab = 219.398,ac =0.17679,ad =  0.00241144;
	float ba =  -56.6164,bbb = 162.647,bc = 0.679974,bd  = 0.00423267,be  = -0.324786;
	float ff1a = -21.878,ff1b = 1146.3, ff1c =  0.0267898, ff1d  =  0.000712903;
	float ff2a = -85.8291,ff2b = 2466.59,ff2c = 0.55722, ff2d  =  0.00159572;
	float ArEhito10=0, ArEhite10=0, ArEhito20=0, ArEhite20m=0, ArEhite20p=0;
	int   ArNhito10=0, ArNhite10=0, ArNhito20=0, ArNhite20m=0, ArNhite20p=0;
	float diff=1.0;
	float x = 0, y= 0, f1 = 0, f2 = 0;
	float _costheta = 0;
	float Ethetacorr = 1;
	float Ephicorr = 1;
	float EMC = inputCluster.getEnergy(); 
	TVector3 CluPos = TVector3(inputCluster.getPosition().x,inputCluster.getPosition().y,inputCluster.getPosition().z); 
	float CluTheta = CluPos.Theta();
	float CluPhi = CluPos.Phi()*5.72957795130823229e+01;	
	
	int NCluHits = inputCluster.hits_size();
	for(int s1 = 0; s1 < NCluHits; s1++)
	{
		auto a_hit = inputCluster.getHits(s1);
		auto cellid=a_hit.getCellID();
		int NLayer = m_decoder->get(cellid, "layer");
		if(NLayer > 20){
			if (NLayer%2 ==0){
				ArNhite10 ++;
				ArEhite10 +=a_hit.getEnergy();
			}
			else if (NLayer%2 ==1){
				ArNhito10 ++;
				ArEhito10 +=a_hit.getEnergy();
			}

		}
		else if(NLayer <= 20){
			if (NLayer%2 ==0){
				if(NLayer ==0){
					ArNhite20m ++;
					ArEhite20m +=a_hit.getEnergy();
				}
				if(NLayer >0){
					ArNhite20p ++;
					ArEhite20p +=a_hit.getEnergy();
				}
			}
			else if (NLayer%2 ==1){
				ArNhito20 ++;
				ArEhito20 +=a_hit.getEnergy();
			}
		}
	}
	while(diff>0.01 && EMC > 0)
	{
		float temp_a=sqrt(EMC*EMC);
		x=exp((-EMC+aaa)/ab)+ac/sqrt(EMC)+ad*EMC;
		y=1/(exp((-EMC+ba)/bbb)+bc/sqrt(EMC)+bd*EMC)+be;
		f1=1/(exp((-EMC+ff1a)/ff1b)+ff1c/sqrt(EMC)+ff1d*EMC);
		f2=1/(exp((-EMC+ff2a)/ff2b)+ff2c/sqrt(EMC)+ff2d*EMC);
		EMC=x*ArEhito20+f1*ArEhite20p+y*ArEhito10+f2*ArEhite10;
		float temp_b=sqrt(EMC*EMC);
		diff=fabs(temp_a-temp_b);
	}
	_costheta=cos(CluTheta);

	if(EMC/inputCluster.getEnergy() < 0.5) EMC = inputCluster.getEnergy();

	Ethetacorr=(50/(0.294596*fabs(_costheta)*fabs(_costheta)-1.58336*fabs(_costheta)+49.9219));

	float ModuleCluPhi = 0;
	if(CluPhi > 17.5)
	{
		ModuleCluPhi = CluPhi - int((CluPhi -17.5)/45)*45;
	}
	else
	{
		ModuleCluPhi = CluPhi - int((CluPhi -17.5)/45)*45 + 45;
	}

	if(ModuleCluPhi>=17.5 && ModuleCluPhi<20.5 ) Ephicorr=(50/(-0.122*ModuleCluPhi*ModuleCluPhi+2.76*ModuleCluPhi+38.04));
	if(ModuleCluPhi>=20.5 && ModuleCluPhi<22.5 ) Ephicorr=(50/(0.22*ModuleCluPhi*ModuleCluPhi-6.43*ModuleCluPhi+81.69));
	if(ModuleCluPhi>=22.5 && ModuleCluPhi<62.5 ) Ephicorr=50*(0.00839675/sqrt(ModuleCluPhi)+0.000369287*logf(ModuleCluPhi)+0.0174033);

	EMC=EMC*Ethetacorr*Ephicorr;

	return EMC;
}

std::vector<float> ArborToolLCIO::ClusterTime(edm4hep::ConstCluster inputCluster)
{
	std::vector<float> CluTimeVector; 
	CluTimeVector.clear();

	int NHit = inputCluster.hits_size();
	float currhitTime = 0;
	float currhitoriTime = 0; 
	TVector3 hitpos; 

	std::vector<float> Time;
	Time.clear(); 
	std::map<int, float> CluHitToTime;
	CluHitToTime.clear();
	std::vector<float> OriginalTime;
	OriginalTime.clear();

	int N_PropTimeHit = 0;

	for(int i = 0; i < NHit; i++)
	{
		auto a_hit =inputCluster.getHits(i);
		hitpos = TVector3(a_hit.getPosition().x,a_hit.getPosition().y,a_hit.getPosition().z);
		if(a_hit.getTime() == 0)
		{
			currhitTime = 1001;	
			currhitoriTime = 1001; 
		}
		else
		{
			currhitTime = a_hit.getTime();
			currhitoriTime = a_hit.getTime();
			N_PropTimeHit++;
		}
		Time.push_back( currhitTime );
		OriginalTime.push_back( currhitoriTime );
		CluHitToTime[i] = currhitTime;			
	}

	std::vector<int> TimeOrder = SortMeasure(Time, 0);
	std::vector<int> OriTimeOrder = SortMeasure(OriginalTime, 0);
	float PeakTime = 0;
	float AverageTime = 0; 
	int NCount = 0; 
	float ptime = 0; 
	int Break = 0; 
	TVector3 StartP, EndP, hitP; 

	for(int j0= 0; j0 < N_PropTimeHit; j0++)
	{
		auto a_hit = inputCluster.getHits(OriTimeOrder[j0]);
		hitP = TVector3(a_hit.getPosition().x,a_hit.getPosition().y,a_hit.getPosition().z);
		if(N_PropTimeHit > 10)
		{
			if( j0 < 5)
				StartP += 0.2*hitP;
			else if(j0 > 4 && j0 < 10)
				EndP += 0.2*hitP;
		}
		else
		{
			if( j0 < N_PropTimeHit/2.)
				StartP += N_PropTimeHit/2.*hitP;
			else
				EndP += N_PropTimeHit/2.*hitP;
		}
	}

	for(int j = 0; j < NHit; j++)
	{
		ptime = CluHitToTime[TimeOrder[j]];

		if(j == 0)
		{
			CluTimeVector.push_back(ptime);
		}
		if(j == NHit - 1)
			CluTimeVector.push_back(ptime);

		if(j > 0)
		{
			//cout<<"ptime "<<ptime<<" cluhittotime "<<CluHitToTime[TimeOrder[j - 1]]<<endl;
			if( (ptime - CluHitToTime[TimeOrder[j - 1]]) < 3)
			{
				if( Break == 0 )
				{
					PeakTime += ptime; 
					NCount ++; 
				}
			}	
			else
			{
				Break = 1;
			}
		}

		AverageTime += ptime; 
	}
	if(NCount)
		PeakTime = float(PeakTime)/NCount; 

	CluTimeVector.push_back(float(AverageTime)/NHit);	
	CluTimeVector.push_back(PeakTime);
	CluTimeVector.push_back(NCount);
	if( N_PropTimeHit > 1 )	// Direction: Cos Theta Between Time Flow And Position
		CluTimeVector.push_back((EndP + StartP).Angle(EndP - StartP));
	else
		CluTimeVector.push_back(1001);

	return CluTimeVector; 
}




