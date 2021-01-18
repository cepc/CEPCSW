#include "BushConnect.hh"
#include "ArborTool.h"
#include "ArborToolLCIO.hh"
#include "DetectorPos.hh"

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
#include "edm4hep/ReconstructedParticle.h"
#include "edm4hep/ClusterCollection.h"
#include "edm4hep/SimCalorimeterHitCollection.h"
#include "edm4hep/MCRecoCaloAssociationCollection.h"
#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/ReconstructedParticleCollection.h"

#include "cellIDDecoder.h"

#include "DD4hep/Detector.h"
#include "DD4hep/IDDescriptor.h"
#include "DD4hep/Plugins.h"

#include <DDRec/DetectorData.h>
#include <DDRec/CellIDPositionConverter.h>
#include "DetInterface/IGeomSvc.h"


#include <values.h>
#include <string>
#include <iostream>
#include <cmath>
#include <stdexcept>
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <Rtypes.h>
#include <sstream>
#include <set>
#include <TVector3.h>
#include <vector>
#include <algorithm>

#include "HelixClassD.hh"		

using namespace std;

//const string ECALCellIDDecoder  = "M:3,S-1:3,I:9,J:9,K-1:6";

DECLARE_COMPONENT(BushConnect)

BushConnect::BushConnect(const std::string& name, ISvcLocator* svcLoc)
	: GaudiAlgorithm(name, svcLoc)
{

}


StatusCode BushConnect::initialize() {

      ISvcLocator* svcloc = serviceLocator();
      m_ArborToolLCIO=new ArborToolLCIO("arborTools",svcloc);
//printParameters();
	//Cluflag.setBit(LCIO::CHBIT_LONG);
	return GaudiAlgorithm::initialize();
}

void BushConnect::Clean(){
	MCPTrack_Type.clear();
	Track_Energy.clear();
	Track_Type.clear();
	Track_Theta.clear();
	Track_EndPoint.clear();
	Track_P3.clear();
	TrackStartPoint.clear();
	SortedTracks.clear();
	ChCoreID.clear();
	chargedclustercore.clear();
	selfmergedcluster.clear();
	non_chargedclustercore.clear();		//all clusters besides charged cluster core
	ClusterType_1stID.clear();
	CluFD.clear();
	CluT0.clear();
	CluCoG.clear();
	Clu_Depth.clear();
}

void BushConnect::TrackSort() //, &std::map<Track*, int>Track_Tpye, &std::map<Track*, float> Track_Energy)
{
	try{
		auto TrackColl = m_trkcol.get();


		int NTrk = TrackColl->size();
		cout<<NTrk<<endl;
		float D0 = 0;
		float Z0 = 0;
		int NTrkHit = 0;
		const float mass = 0.139;	//Pion Mass
		TVector3 EndPointPos, StartPointPos; 
		int TrackType = 0; 

		std::vector<edm4hep::ConstTrack> tracks_HQ_Barrel; 
		std::vector<edm4hep::ConstTrack> tracks_HQ_Endcap;
		std::vector<edm4hep::ConstTrack> tracks_HQ_Shoulder;
		std::vector<edm4hep::ConstTrack> tracks_HQ_Forward; 
		std::vector<edm4hep::ConstTrack> tracks_MQ_Barrel;
		std::vector<edm4hep::ConstTrack> tracks_MQ_Endcap;
		std::vector<edm4hep::ConstTrack> tracks_MQ_Shoulder;
		std::vector<edm4hep::ConstTrack> tracks_MQ_Forward;
		std::vector<edm4hep::ConstTrack> tracks_Vtx; 
		std::vector<edm4hep::ConstTrack> tracks_LQ; 
		std::vector<edm4hep::ConstTrack> tracks_LE; 
		std::vector<edm4hep::ConstTrack> curr_tracks;
	
		Track_EndPoint.clear();
	
		tracks_HQ_Barrel.clear();
		tracks_HQ_Endcap.clear();
		tracks_HQ_Shoulder.clear();
		tracks_HQ_Forward.clear();
		tracks_MQ_Barrel.clear();
		tracks_MQ_Endcap.clear();
		tracks_MQ_Shoulder.clear();
		tracks_MQ_Forward.clear();
		tracks_Vtx.clear();
		tracks_LQ.clear();
		tracks_LE.clear();

		std::vector<edm4hep::ConstTrack> tracks_ILL;
		tracks_ILL.clear();
		std::vector<edm4hep::ConstTrack> tracks_preInteraction;
		tracks_preInteraction.clear();	//Used to denote pion and electron interaction inside TPC/Tracker. Simply vetoed for avoid double counting... but muon may still be problematic. Better way of treating would be find the cascade photons & tracks - clusters, and veto all the daughters instead of mother. Similar can done for Kshort...
	// Condition, tracks_head to others tail. head position far from boundary. and, track energy >= sum of cascade

		std::vector<int> TrackOrder; 
		TrackOrder.clear();	
		std::map<edm4hep::ConstTrack, int> Track_Index; 
		Track_Index.clear();
		Track_Energy.clear();
		Track_Type.clear();
		Track_P3.clear();
		Track_EndPoint.clear();
		TrackStartPoint.clear();
	
		for(int i0 = 0; i0 < NTrk; i0++)
		{
			auto a_Trk = (*TrackColl)[i0];
			NTrkHit = a_Trk.getTrackerHits().size();		
			EndPointPos = TVector3((a_Trk.getTrackerHits(NTrkHit - 1)).getPosition().x,(a_Trk.getTrackerHits(NTrkHit - 1)).getPosition().y,(a_Trk.getTrackerHits(NTrkHit - 1)).getPosition().z);	
			StartPointPos = TVector3((a_Trk.getTrackerHits(0)).getPosition().x,(a_Trk.getTrackerHits(0)).getPosition().y,(a_Trk.getTrackerHits(0)).getPosition().z);
			Track_EndPoint[a_Trk] = EndPointPos;
			TrackStartPoint[a_Trk] = StartPointPos;
	
			HelixClassD * TrkInit_Helix = new HelixClassD();
			TrkInit_Helix->Initialize_Canonical(a_Trk.getTrackStates(0).phi, a_Trk.getTrackStates(0).D0, a_Trk.getTrackStates(0).Z0, a_Trk.getTrackStates(0).omega, a_Trk.getTrackStates(0).tanLambda, BField);
			float TrackEn = mass*mass;
	
			for (int q3 = 0; q3 < 3; q3 ++)
			{
				TrackEn += (TrkInit_Helix->getMomentum()[q3])*(TrkInit_Helix->getMomentum()[q3]);
			}
			TVector3 TrkMom(TrkInit_Helix->getMomentum()[0],TrkInit_Helix->getMomentum()[1],TrkInit_Helix->getMomentum()[2]);
	
			TrackEn = sqrt(TrackEn);
			Track_Energy[a_Trk] = TrackEn;
			Track_Theta[a_Trk] = TrkMom.Theta();
			// Track_Phi[a_Trk] = TrkMom.Phi();
			Track_P3[a_Trk] = TrkMom;		
	
			delete TrkInit_Helix;
		}
	
		TVector3 currEp, currSp;
		float currMotherEn = 0;
		float sumDauEn = 0; 
	
		for(int i1 = 0; i1 < NTrk; i1++)
		{
			auto a_Trk = (*TrackColl)[i1];
			currEp = Track_EndPoint[a_Trk];
	
			if( currEp.Perp() < 1600 && currEp.Perp() > 400 && abs(currEp.Z()) < 2000 )	//Only check 
			{
				currMotherEn = Track_Energy[a_Trk];
				sumDauEn = 0;	
				for(int i2 = 0; i2 < NTrk; i2++)
				{
					auto b_Trk = (*TrackColl)[i2];
					if(i2 != i1)
					{
						currSp = TrackStartPoint[b_Trk];
						if( (currEp - currSp).Mag() < 40  )
							sumDauEn += Track_Energy[b_Trk];
					}
				}
				if(currMotherEn + 0.1 > 0.9*sumDauEn && currMotherEn > 3 && sumDauEn > 0 )	//Some protection is always needed...
				{
					tracks_preInteraction.push_back(a_Trk);
				}
			}
		}
	
		for(int t0 = 0; t0 < NTrk; t0++)
		{
			auto a_Trk = (*TrackColl)[t0];
			D0 = a_Trk.getTrackStates(0).D0;
			Z0 = a_Trk.getTrackStates(0).Z0;
			NTrkHit = a_Trk.getTrackerHits().size();
			auto last_hit = a_Trk.getTrackerHits(NTrkHit - 1);
			EndPointPos = TVector3(last_hit.getPosition().x,last_hit.getPosition().y,last_hit.getPosition().z);
			Track_EndPoint[a_Trk] = EndPointPos;
			StartPointPos = TVector3((a_Trk.getTrackerHits(0)).getPosition().x,(a_Trk.getTrackerHits(0)).getPosition().y,(a_Trk.getTrackerHits(0)).getPosition().z);
			Track_Index[a_Trk] = t0;
	
			if( NTrkHit > 9 || (fabs(EndPointPos.Z()) > LStar - 500 && EndPointPos.Perp() < TPCInnerRadius ) || fabs(EndPointPos.Z()) > ECALHalfZ - 200  )		// Min requirement for track quality
			{	// LStar - 500, suppose to be the last Disk Position
	
				if( find(tracks_preInteraction.begin(), tracks_preInteraction.end(), a_Trk ) != tracks_preInteraction.end() )
				{
					cout<<"So We Drop it! "<<Track_Energy[a_Trk]<<endl; 
					continue; 
				}
	
				TrackType = 0;
	
				if((Track_Energy[a_Trk] < 1.0 && fabs(Track_Theta[a_Trk]-1.57)< 0.4) || (fabs(Track_Theta[a_Trk]-1.57) >= 0.4 && log10(Track_Energy[a_Trk]) < -(fabs(Track_Theta[a_Trk]-1.57)-0.4)*0.2/0.3 ))
				{
					TrackType = 100;
				}
				else if( fabs(EndPointPos.Z()) > ECALHalfZ - 500 && EndPointPos.Perp() > TPCOuterRadius - 300  )	//Shoulder
				{
					TrackType = 30;
				}
				else if( fabs(EndPointPos.Z()) > LStar - 500 && EndPointPos.Perp() < TPCInnerRadius )		//Forward
				{
					TrackType = 40;
				}
				else if( EndPointPos.Perp() > TPCOuterRadius - 100 )		//Barrel
				{
					TrackType = 10;
				}
				else if( fabs(EndPointPos.Z()) > ECALHalfZ - 200 )		//Endcap
				{
					TrackType = 20; 
				}
	
				if( fabs(D0) < 1 && fabs(Z0) < 1 )
				{
					TrackType += 1;
				}
	
				Track_Type[a_Trk] = TrackType; 
	
				if(TrackType == 11)
					tracks_HQ_Barrel.push_back(a_Trk);
				else if(TrackType == 21)
					tracks_HQ_Endcap.push_back(a_Trk);
				else if(TrackType == 31)
					tracks_HQ_Shoulder.push_back(a_Trk);
				else if(TrackType == 41)
					tracks_HQ_Forward.push_back(a_Trk);
				else if(TrackType == 10)
					tracks_MQ_Barrel.push_back(a_Trk);
				else if(TrackType == 20)
					tracks_MQ_Endcap.push_back(a_Trk);
				else if(TrackType == 30)
					tracks_MQ_Shoulder.push_back(a_Trk);
				else if(TrackType == 40)
					tracks_MQ_Forward.push_back(a_Trk);
				else if(TrackType == 1)
					tracks_Vtx.push_back(a_Trk);
				else if(TrackType == 101)
					tracks_LE.push_back(a_Trk);
				else if( (StartPointPos.Mag() > 50 && EndPointPos.Mag() < 1000 && NTrkHit < 50) || TrackType == 100  )
				{			tracks_ILL.push_back(a_Trk);
				}
				else
					tracks_LQ.push_back(a_Trk);
			}
		}
		cout<<"LQ"<<tracks_LQ.size()<<" "<<tracks_ILL.size()<<endl;
	
		std::vector<float > currTrkMomentum;
		std::vector<int> currTrkIndex;
	
		for(int t1 = 0; t1 < 11; t1++)
		{
			currTrkMomentum.clear();
			currTrkIndex.clear();
			curr_tracks.clear();
			if(t1 == 0)
				curr_tracks = tracks_HQ_Endcap;
			else if(t1 == 1)
				curr_tracks = tracks_HQ_Barrel;
			else if(t1 == 2)
				curr_tracks = tracks_MQ_Endcap;
			else if(t1 == 3)
				curr_tracks = tracks_MQ_Barrel;
			else if(t1 == 4)
				curr_tracks = tracks_HQ_Shoulder;
			else if(t1 == 5)
				curr_tracks = tracks_MQ_Shoulder;
			else if(t1 == 6)
				curr_tracks = tracks_HQ_Forward;
			else if(t1 == 7)
				curr_tracks = tracks_MQ_Forward;
			else if(t1 == 8)
				curr_tracks = tracks_Vtx;
			else if(t1 == 9)			
				curr_tracks = tracks_LQ; 
			else if(t1 == 10)			
				curr_tracks = tracks_LE; 
	
	
			int N_currTrack = curr_tracks.size();
	
			for(int t2 = 0; t2 < N_currTrack; t2++)
			{
				auto tmpTrk = curr_tracks[t2];
				currTrkMomentum.push_back(Track_Energy[tmpTrk]);
			}
	
			currTrkIndex = SortMeasure(currTrkMomentum, 1);
	
			for(int t3 = 0; t3 < N_currTrack; t3++)
			{
				auto b_tmpTrk = curr_tracks[currTrkIndex[t3]];
				if(t1 < 9 || Track_Energy[b_tmpTrk] < 10)
					TrackOrder.push_back(Track_Index[b_tmpTrk]);
			}
		}
	
		for(unsigned int t4 = 0; t4 < TrackOrder.size(); t4++)
		{
			auto b_Trk =  (*TrackColl)[t4];
			SortedTracks.push_back(b_Trk);
		}
	}catch(GaudiException &e){}
}

void BushConnect::BushSelfMerge()
{
	auto CaloClu = m_clucol.get();
	int NClu = CaloClu->size();

	std::vector<edm4hep::ConstCluster > Core_1st; 
	std::vector<edm4hep::ConstCluster > Frag_1st;
	std::vector<edm4hep::ConstCluster > UnId_1st; 
	Core_1st.clear();
	Frag_1st.clear();
	UnId_1st.clear();

	float CluDepth = 0; 
	float CluEn = 0;
	int CluSize = 0; 
	TVector3 PosCluSeed, PosSeedDiff, PosCoGDiff, PosSeedA, PosSeedB; 

	int NJoints = 0; 	
	int SmallCluSize = 0; float DeeperDepth = 0; float LaterT0 = 0; 
	float Depth_A = 0; float Depth_B = 0;
	int Size_A = 0; int Size_B = 0; 

	TMatrixF FlagMerge(NClu, NClu);

	cout<<NClu<<" clusters"<<endl;
	for(int i0 = 0; i0 < NClu; i0++)
	{
		auto a_clu = (*CaloClu)[i0];
		CluFD[a_clu] = m_ArborToolLCIO->FDV3(a_clu);
		CluT0[a_clu] = m_ArborToolLCIO->ClusterT0(a_clu);
		CluCoG[a_clu] = m_ArborToolLCIO->ClusterCoG(a_clu);
		PosCluSeed = TVector3(a_clu.getPosition().x,a_clu.getPosition().y,a_clu.getPosition().z);
		Clu_Depth[a_clu] = DisSeedSurface(PosCluSeed);
	}
	// CluCoG_Top20Percent Might be used to improve Photon Split Remerge performance

	for(int s0 = 0; s0 < NClu; s0++)
	{
		auto a_clu = (*CaloClu)[s0];
		PosSeedA = TVector3(a_clu.getPosition().x,a_clu.getPosition().y,a_clu.getPosition().z);
		Depth_A = Clu_Depth[a_clu];
		Size_A = a_clu.hits_size();

		for(int s1 = s0 + 1; s1 < NClu; s1++)
		{
			auto b_clu = (*CaloClu)[s1];
			NJoints = m_ArborToolLCIO->JointsBetweenBush(a_clu, b_clu, 10);
			PosSeedB = TVector3(b_clu.getPosition().x,b_clu.getPosition().y,b_clu.getPosition().z);
			Depth_B = Clu_Depth[b_clu];
			Size_B = b_clu.hits_size();

			PosSeedDiff = PosSeedA - PosSeedB;
			DeeperDepth = std::max(Depth_A, Depth_B);
			LaterT0 = std::max(CluT0[a_clu], CluT0[b_clu]);
			PosCoGDiff = CluCoG[a_clu] - CluCoG[b_clu];

			if(NJoints&& PosSeedDiff.Mag() < 40 + 0.05*DeeperDepth )	//And depth...
			{
				SmallCluSize = std::min( Size_A, Size_B );

				//if( SmallCluSize < 10 || LaterT0 > 10 || DeeperDepth > 40 || (NJoints > 8 && PosCoGDiff.Mag()*PosCoGDiff.Angle(PosSeedB) < 15 ))
				if( SmallCluSize < 10 || LaterT0 > 10 || DeeperDepth > 40 || (NJoints > 8 && PosCoGDiff.Mag()*PosCoGDiff.Angle(PosSeedB) < 15 ) || NJoints>16)

					//if( SmallCluSize < 10 || LaterT0 > 10 || DeeperDepth > 40 || (NJoints > 8 ) )
				{	
					FlagMerge[s0][s1] = 1.0;
					FlagMerge[s1][s0] = 1.0;
				}
			}
			//Head Tail Connection. Could be more sophsticate && should be very strict.
			if( PosSeedA.Angle(PosSeedB) < 0.1 && PosSeedDiff.Mag() < 1000 && PosSeedDiff.Mag()*PosSeedA.Angle(PosSeedB) < 60 + 0.02*DeeperDepth && ((CluFD[a_clu] < 0.2 && Size_A > 6) || (CluFD[b_clu] < 0.2 && Size_B > 6)) )
			{
				if( (PosSeedA.Mag() > PosSeedB.Mag() && PosSeedA.Angle(PosSeedB - PosSeedA) < 0.2) || (PosSeedB.Mag() > PosSeedA.Mag() && PosSeedA.Angle(PosSeedA - PosSeedB) < 0.2) )
				{
					FlagMerge[s0][s1] = 2.0;
					FlagMerge[s1][s0] = 2.0;
				}
			}
		}
	}

	std::vector<edm4hep::ConstCluster> OriInputEHBushes = m_ArborToolLCIO->CollClusterVec(CaloClu);
	TMatrixF MergeSYM = MatrixSummarize(FlagMerge);
	auto CloseMergedCaloClu = m_ArborToolLCIO->ClusterVecMerge( OriInputEHBushes, MergeSYM, clucol_merged);

	std::map<edm4hep::ConstCluster,float> MinDisSeedToBush;
	MinDisSeedToBush.clear();
	for(int i0 = 0; i0 < CloseMergedCaloClu->size(); i0++)
	{
		auto a_clu = (*CloseMergedCaloClu)[i0];
		PosCluSeed = TVector3(a_clu.getPosition().x,a_clu.getPosition().y,a_clu.getPosition().z);
		float tmpmindis = 1e10;
		for(int i1 = 0; i1 < CloseMergedCaloClu->size(); i1++)
		{
			auto b_clu = (*CloseMergedCaloClu)[i1];
			if(i1 != i0)
			{
				if(m_ArborToolLCIO->DisPointToBush(PosCluSeed,b_clu) < tmpmindis) 
					tmpmindis = m_ArborToolLCIO->DisPointToBush(PosCluSeed,b_clu);  
			}
		}
		MinDisSeedToBush[a_clu] = tmpmindis;
	}

	for(int i0 = 0; i0 < CloseMergedCaloClu->size(); i0++)
	{
		auto a_clu = (*CloseMergedCaloClu)[i0];
		PosCluSeed = TVector3(a_clu.getPosition().x,a_clu.getPosition().y,a_clu.getPosition().z);
		CluDepth = DisSeedSurface(PosCluSeed);
		CluEn = a_clu.getEnergy();
		CluSize = a_clu.hits_size();

		if( CluSize > 10 && ( (CluEn > 2.0 + 0.002*CluDepth && CluDepth < 22 ) || CluEn > 5.0) )
		{
			Core_1st.push_back(a_clu);
		}
		else if( (CluSize > 10 || CluEn > 0.2) && MinDisSeedToBush[a_clu] > 20 && m_ArborToolLCIO->ClusterT0(a_clu) < 1.2 ) // && (PhotonTag(a_clu) == 1 || ClusterFlag1st(a_clu) == 11  ))
		{
			Core_1st.push_back(a_clu);
		}
		else if( CluSize < 5 && CluEn < 0.3 && CluDepth > 40 )
		{
			Frag_1st.push_back(a_clu);
		}
		else
		{
			UnId_1st.push_back(a_clu);
		}
	}

	std::vector<edm4hep::ConstCluster > UndefFrag_1stAB = m_ArborToolLCIO->ClusterAbsorbtion(UnId_1st, Frag_1st, 50, 0.02);
	selfmergedcluster = m_ArborToolLCIO->ClusterAbsorbtion(Core_1st, UndefFrag_1stAB, 50, 0.02);
	auto CluAB_1st=m_ArborToolLCIO->ClusterVecColl(selfmergedcluster,m_1stclucol);
}

void BushConnect::TagCore() 
{
	int NTrk = SortedTracks.size();
	int NClu = selfmergedcluster.size();
	int currTrackType = 0;
	float currTrkEn = 0;
	float DisMatrix_Track_Clu_E[NTrk][NClu];
	float TimeMatrix_Track_Clu_E[NTrk][NClu];
	float CluDepth = 0;
	float CoreMergeDistanceDepthCorrector = 0;

	TVector3 TrkEndPoint(0, 0, 0);	
	TVector3 CluPos(0, 0, 0);
	std::map<edm4hep::ConstCluster, int> BushTouchFlag; 
	std::map<edm4hep::ConstTrack, edm4hep::ConstCluster> FCMap_Track_CHCore;
	std::map<edm4hep::ConstTrack, std::vector<edm4hep::ConstCluster>> FCMap_Track_CHCore_new;
	std::map<int, int> Closest_Trk_Clu_Map;
	std::vector<edm4hep::ConstCluster> TightLinkedCluster;
	TightLinkedCluster.clear();
	Closest_Trk_Clu_Map.clear();
	BushTouchFlag.clear();
	FCMap_Track_CHCore.clear();
	FCMap_Track_CHCore_new.clear();

	Clu_Depth.clear();

	for(int s0 = 0; s0 < NTrk; s0++)
	{
		for(int s1 = 0; s1 < NClu; s1++)
		{
			DisMatrix_Track_Clu_E[s0][s1] = 1.0E10;
			TimeMatrix_Track_Clu_E[s0][s1] = 1.0E10; 
		}
	}

	for(int t0 = 0; t0 < NClu; t0++)
	{
		auto a_clu = selfmergedcluster[t0];
		TVector3 PosCluSeed = TVector3(a_clu.getPosition().x,a_clu.getPosition().y,a_clu.getPosition().z);
		Clu_Depth[a_clu] = DisSeedSurface(PosCluSeed);
	}

	//~~~~~~~ find the closest cluster first...

	for(int g0 = 0; g0 < NTrk; g0++)
	{
		auto a_trk = SortedTracks[g0];
		float ClosestDis = 1.0E9;  
		int ClosestCluIndex = -1; 
		int ClosestNC = 1E9;
		float TrackEn = Track_Energy[a_trk];
		TrkEndPoint = Track_EndPoint[a_trk];

		currTrackType = Track_Type[a_trk];
		TVector3 TrkP3 = Track_P3[a_trk];

		for(int g1 = 0; g1 < NClu; g1++)
		{
			auto fccand_bush = selfmergedcluster[g1];
			float* Dis = m_ArborToolLCIO->SimpleDisTrackClu(a_trk, fccand_bush);
			float Time = m_ArborToolLCIO->SimpleBushTimeTrackClu(a_trk, fccand_bush);
			int NC = m_ArborToolLCIO->SimpleBushNC(a_trk, fccand_bush);
			if(Dis[2] > -0.1)
			{
				DisMatrix_Track_Clu_E[g0][g1] = Dis[2];
				TimeMatrix_Track_Clu_E[g0][g1] = Time;
				if( Dis[2] < ClosestDis ) // && ThetaDiff < 0.05 + 0.1/Track_Energy[a_trk] )
				{
					ClosestDis = Dis[2]; 
					ClosestCluIndex = g1;
					ClosestNC = NC;
				}
			}
		}

		//Diag for mimic
		cout<<" Z R "<<TrkEndPoint.Z()<<" MM "<<TrkEndPoint.Perp()<< " ClosestCluIndex "<<ClosestCluIndex<<" ClosestDis "<<ClosestDis <<endl; 
		//End Diag

		if( ClosestDis < 15 + 15./TrackEn && ClosestCluIndex > -0.1 && (ClosestNC < 3 || abs(TrkP3.Theta() - 1.57) < 0.01 ) ) 
		{
			auto candiclu= selfmergedcluster[ClosestCluIndex];
			CluPos = TVector3(candiclu.getPosition().x,candiclu.getPosition().y,candiclu.getPosition().z);
			float TrackEndPDis = (TrkEndPoint - CluPos).Mag();
			if( (TrackEndPDis < 400 + 400./TrackEn || (TrackEn < 1 && candiclu.getEnergy() < 2 )) && ( fabs(TrackEn - candiclu.getEnergy() ) < 3.0*sqrt(TrackEn) + 1.0 ||  candiclu.getEnergy() < 3) && (TrkEndPoint.Z()*CluPos.Z() > 0 || abs(TrkEndPoint.Z() - CluPos.Z()) < 100 ) ) // && AngDiff < 0.4 )
			{
				Closest_Trk_Clu_Map[g0] = ClosestCluIndex;
				BushTouchFlag[candiclu] = 0;
			}
		}
	}	//~~~~~~~ end of finding closest cluster

	cout<<"NTrk "<<NTrk<<endl;
	for(int i0 = 0; i0 < NTrk; i0++)  //Dropped Size can exist
	{
		auto a_trk = SortedTracks[i0];
		currTrackType = Track_Type[a_trk];
		currTrkEn = Track_Energy[a_trk];

		TrkEndPoint = Track_EndPoint[a_trk];
		FCMap_Track_CHCore_new[a_trk].clear();
		if( Closest_Trk_Clu_Map.find(i0) != Closest_Trk_Clu_Map.end() )
		{
			auto closeClu = selfmergedcluster[Closest_Trk_Clu_Map[i0]];
			FCMap_Track_CHCore_new[a_trk].push_back(closeClu);
		}

		for(int j0 = 0; j0 < NClu; j0++)
		{
			auto fccand_bush = selfmergedcluster[j0]; 
			float Dis = DisMatrix_Track_Clu_E[i0][j0]; //SimpleDisTrackClu(a_trk, fccand_bush);
			float BushTime = TimeMatrix_Track_Clu_E[i0][j0];
			CluDepth = Clu_Depth[fccand_bush];
			CluPos = TVector3(fccand_bush.getPosition().x,fccand_bush.getPosition().y,fccand_bush.getPosition().z);
			int currCluType = m_ArborToolLCIO->ClusterFlag1st(fccand_bush);

			CoreMergeDistanceDepthCorrector = 0;
			if(CluDepth > 20)
				CoreMergeDistanceDepthCorrector = 20;
			else if(CluDepth > 10)
				CoreMergeDistanceDepthCorrector = 10;

			float ProjectiveDis_TrackEndP_CluPos = (CluPos - TrkEndPoint).Mag()*(CluPos - TrkEndPoint).Angle(CluPos);
			if(log10(BushTime) < 3.5 && BushTime > 0 && currTrackType != 101 && Dis < 10 + CoreMergeDistanceDepthCorrector && Dis > -0.1 && BushTouchFlag.find(fccand_bush) == BushTouchFlag.end() && (currTrkEn > 3 || fccand_bush.getEnergy() < 5 || currCluType == 13 ) && ProjectiveDis_TrackEndP_CluPos < 100 ) // && (TrackEndPDis < 400 || currTrackType != 11) && (CluPos.Z()*TrkEndPoint.Z() > 0 || fabs((CluPos-TrkEndPoint).Z()) < 20 ) && (fccanden - currTrkEn <  something)  )
			{
				FCMap_Track_CHCore_new[a_trk].push_back(fccand_bush);
				BushTouchFlag[fccand_bush] = currTrackType;

			}
		}
		//Diag 4 MIMIC
	}

	//Added by Dan for multi track linked to one cluster
	for(int i0 = 0; i0 < NTrk; i0++)  //Dropped Size can exist
	{
		auto a_trk = SortedTracks[i0];
		TrkEndPoint = Track_EndPoint[a_trk];
		int naclu=FCMap_Track_CHCore_new[a_trk].size();


		for(int i1=0;i1<naclu;i1++)
		{
			auto aclu=FCMap_Track_CHCore_new[a_trk][i1];
			float dis_a=DisMatrix_Track_Clu_E[i0][i1];
			for(int j0=i0+1; j0<NTrk; j0++)
			{
				int breakflag=0;
				auto b_trk = SortedTracks[j0];
				int nbclu=FCMap_Track_CHCore_new[b_trk].size();

				for(int j1=0;j1<nbclu;j1++)
				{
					auto bclu=FCMap_Track_CHCore_new[b_trk][j1];
					float dis_b=DisMatrix_Track_Clu_E[j0][j1];
					if(aclu==bclu){
						if(dis_a>=dis_b){
							FCMap_Track_CHCore_new[a_trk].erase(FCMap_Track_CHCore_new[a_trk].begin()+i1);
							breakflag=1;
							break;
						}
						else{
							FCMap_Track_CHCore_new[b_trk].erase(FCMap_Track_CHCore_new[b_trk].begin()+j1);
						}
					}
				}
				if(breakflag)break;
			}
		}
		if(FCMap_Track_CHCore_new[a_trk].size() > 0 ) // && EcalCoreEnergy + HcalCoreEnergy < 2.0*currTrkEn )...
		{
			auto chcorecluster_eh =  m_ArborToolLCIO->NaiveMergeClu(FCMap_Track_CHCore_new[a_trk]);
			edm4hep::ConstCluster chcorecluster_ehCon=chcorecluster_eh;
			FCMap_Track_CHCore[a_trk] = chcorecluster_ehCon;
			chargedclustercore.push_back(chcorecluster_ehCon);
		}
	}// End by Dan

	// Diagnosis in Tagcore
	for(int t0 = 0; t0 < NClu; t0++)
	{
		auto a_clu = selfmergedcluster[t0];
		if( BushTouchFlag.find(a_clu) == BushTouchFlag.end() ) //Might be a neutral core
		{
			non_chargedclustercore.push_back(a_clu);
		}
	}

	int Track_Core_ID = -99;

	edm4hep::ReconstructedParticleCollection* chargeparticleCol = m_chargeparticleCol.createAndPut();
	edm4hep::ClusterCollection* chargedcoreclusterCol = m_chargedcoreclusterCol.createAndPut();
	for(int j5 = 0; j5 < NTrk; j5++)
	{
		auto a_trk = SortedTracks[j5];
		Track_Core_ID = -99;
		currTrackType = Track_Type[a_trk];
		auto chargeparticle=chargeparticleCol->create();
		chargeparticle.setEnergy( Track_Energy[a_trk] );
		chargeparticle.setCharge(a_trk.getTrackStates(0).omega/fabs(a_trk.getTrackStates(0).omega));
		TVector3 Ptrack = Track_P3[a_trk];
		edm4hep::Vector3f currTrkP = edm4hep::Vector3f( Ptrack.X(), Ptrack.Y(), Ptrack.Z() );
		chargeparticle.setMomentum(currTrkP); 
		chargeparticle.addToTracks( a_trk );
		auto a_clu = FCMap_Track_CHCore[a_trk];
		if( FCMap_Track_CHCore[a_trk].hits_size()>0 )		// No really need to pertect, as quality will be controled in Part.Reco
		{
			edm4hep::Cluster chargedcorecluster = chargedcoreclusterCol->create();

			m_ArborToolLCIO->NaiveCluConst(FCMap_Track_CHCore[a_trk],chargedcorecluster);
			edm4hep::ConstCluster chargedcoreclusterCon=chargedcorecluster;
			chargeparticle.addToClusters(chargedcoreclusterCon);
			Track_Core_ID = m_ArborToolLCIO->ClusterFlag(a_clu, a_trk);
		}
		chargeparticle.setType(Track_Core_ID);
		ChCoreID[chargeparticle] = Track_Core_ID;
	}


}

void BushConnect::ParticleReco()
{

	auto col_IsoHit = m_col_IsoHit.get();
	std::vector<edm4hep::ConstCalorimeterHit> IsoHits = m_ArborToolLCIO->CollHitVec(col_IsoHit, 0);

	edm4hep::ReconstructedParticleCollection* arborrecoparticleCol = m_arborrecoparticleCol.createAndPut();

	edm4hep::ClusterCollection* mergedclu_chCol = m_mergedclu_chCol.createAndPut();
	edm4hep::ClusterCollection* mergedclu_neCol = m_mergedclu_neCol.createAndPut();

	auto ChargedCore = m_chargeparticleCol.get(); 
	int NChargedObj = ChargedCore->size();
	int NNeutralCluster = non_chargedclustercore.size();
	double DisMatrix_Core_Neutral[NChargedObj][NNeutralCluster][2];		//Define different types of distances; 

	float CluDepth = 0;
	std::map<edm4hep::ConstCluster, double> CluDepthMap; 
	CluDepthMap.clear();
	int currChargeCoreType = 0;  
	TVector3 CluPos; 

	std::vector<edm4hep::ConstCluster> loosecandicluster; 
	std::vector<edm4hep::ConstCluster> tightcandicluster;		//Muon potential candi?
	std::vector<edm4hep::ConstCluster> mergedcluster; 			//tmp for each charged P
	std::vector<edm4hep::ConstCluster> chargedclustercore_merged; 	//overall
	chargedclustercore_merged.clear();

	std::vector<double> reftightdis; 
	std::vector<double> refloosedis; 
	std::map<edm4hep::ConstCluster, int> NNCTouchFlag; 
	std::vector<edm4hep::ConstTrack> SecondIterTracks;
	SecondIterTracks.clear();

	TVector3 currTrkEnd, neighbourTrkEnd, LeadP; 

	for(int i = 0; i < NChargedObj; i++)
	{
		auto a_recoP_ch = (*ChargedCore)[i];

		loosecandicluster.clear();
		tightcandicluster.clear();
		mergedcluster.clear();
		reftightdis.clear();
		refloosedis.clear();
		auto a_chargedTrk = a_recoP_ch.getTracks(0);
		currTrkEnd = Track_EndPoint[a_chargedTrk];
		currChargeCoreType = ChCoreID[a_recoP_ch];
		int currTrkType = Track_Type[a_chargedTrk];

		float CurrClusterEnergy = 0;
		float CurrTrackEnergy = Track_Energy[a_chargedTrk];
		edm4hep::ConstCluster a_chargedClu;
		if(a_recoP_ch.clusters_size() != 0)
		{
			a_chargedClu = a_recoP_ch.getClusters(0);
			CurrClusterEnergy = a_chargedClu.getEnergy();
			mergedcluster.push_back(a_chargedClu);		//Actually can use this chance to question if previous energy are balance...
		}

		float MinDisToNoClusterTrk = 1.0E10; 
		float MinDisToOtherTrack = 1.0E10;

		for( int is = 0; is < NChargedObj; is++ )
		{
			if(is != i)
			{
				auto b_recoP_ch = (*ChargedCore)[is];
				auto a_neighbourTrk = b_recoP_ch.getTracks(0);
				neighbourTrkEnd = Track_EndPoint[a_neighbourTrk];
				float currDD = (neighbourTrkEnd - currTrkEnd).Mag();
				if( currDD < MinDisToOtherTrack )
				{
					MinDisToOtherTrack = currDD;
				}
			}
		}

		for(int j = 0; j < NNeutralCluster; j++)
		{
			auto a_NeCandiClu = non_chargedclustercore[j];
			float NeCandEn = a_NeCandiClu.getEnergy(); 
			CluPos = TVector3(a_NeCandiClu.getPosition().x,a_NeCandiClu.getPosition().y,a_NeCandiClu.getPosition().z);
			CluDepth = DisSeedSurface(CluPos);
			CluDepthMap[a_NeCandiClu] = CluDepth; 	

			if( ClusterType_1stID[a_NeCandiClu] == 22)   continue; 

			for(int k = 0; k < 2; k++)
			{
				DisMatrix_Core_Neutral[i][j][k] = 1.0E9;
			}

			if(CurrClusterEnergy > 1E-6)	//put by hand...
			{
				DisMatrix_Core_Neutral[i][j][0] = m_ArborToolLCIO->BushDis(a_chargedClu, a_NeCandiClu);
			}
			float* Dis = m_ArborToolLCIO->SimpleDisTrackClu(a_chargedTrk, a_NeCandiClu);
			DisMatrix_Core_Neutral[i][j][1] = Dis[2];

			if( NNCTouchFlag.find(a_NeCandiClu) == NNCTouchFlag.end() && ( currChargeCoreType == 0 || DisMatrix_Core_Neutral[i][j][0] < 1000 ) && currTrkType != 101)
			{			
				if( currChargeCoreType == 130 )			//Matched Muon, should ignore
				{
					if( DisMatrix_Core_Neutral[i][j][1] < 0.2*CluDepth && CluDepth > 200  )	//&& FD?
					{
						tightcandicluster.push_back(a_NeCandiClu);
						reftightdis.push_back( DisMatrix_Core_Neutral[i][j][1] );	//dependence on Cluster Flag & Clu Depth. use some more fancy sort algorithm...
					}
				}
				else if( currChargeCoreType == 131 )
				{
					if( DisMatrix_Core_Neutral[i][j][1] < 0.3*CluDepth && CluDepth > 150 )
					{
						tightcandicluster.push_back(a_NeCandiClu);
						reftightdis.push_back( DisMatrix_Core_Neutral[i][j][1] );	
					}
					else if( DisMatrix_Core_Neutral[i][j][1] < 0.5*CluDepth && CluDepth > 100 )
					{
						loosecandicluster.push_back(a_NeCandiClu);
						refloosedis.push_back( DisMatrix_Core_Neutral[i][j][1] );
					}
				}	
				else if( currChargeCoreType == 110  )		// Electron
				{
					if( DisMatrix_Core_Neutral[i][j][0] < 0.15*CluDepth + 15 )
					{
						tightcandicluster.push_back(a_NeCandiClu);
						reftightdis.push_back( DisMatrix_Core_Neutral[i][j][0] );
					}
				}			
				else if( currChargeCoreType == 111 )		// look behind... might be pion...
				{
					if( DisMatrix_Core_Neutral[i][j][0] < 0.1*CluDepth + 15 && DisMatrix_Core_Neutral[i][j][1] < 0.1*CluDepth + 10 )	//Define Brems Photon region for correct
					{
						tightcandicluster.push_back(a_NeCandiClu);
						if(DisMatrix_Core_Neutral[i][j][0] < DisMatrix_Core_Neutral[i][j][1])	// not fully adequate.
						{
							reftightdis.push_back( DisMatrix_Core_Neutral[i][j][0] );
						}
						else
						{
							reftightdis.push_back( DisMatrix_Core_Neutral[i][j][1] );
						}
					}
					else if( DisMatrix_Core_Neutral[i][j][0] < 0.2*CluDepth + 15 || DisMatrix_Core_Neutral[i][j][1] < 0.2*CluDepth + 15  )
					{	
						loosecandicluster.push_back(a_NeCandiClu);

						if(DisMatrix_Core_Neutral[i][j][0] < DisMatrix_Core_Neutral[i][j][1])   // not fully adequate.
						{
							refloosedis.push_back( DisMatrix_Core_Neutral[i][j][0] );
						}
						else
						{
							refloosedis.push_back( DisMatrix_Core_Neutral[i][j][1] );
						}
					}
				}
				else if( currChargeCoreType == 211 )	//Main Cluster distance oriented
				{
					if(DisMatrix_Core_Neutral[i][j][0] < 0.2*CluDepth)
					{
						tightcandicluster.push_back(a_NeCandiClu);
						reftightdis.push_back( DisMatrix_Core_Neutral[i][j][0] );
					}
				}
				else if( currChargeCoreType == 212 )	//Non_Matched
				{
					if(DisMatrix_Core_Neutral[i][j][0] < 10 + 0.5*CluDepth )	//Energy Dependence...
					{
						tightcandicluster.push_back(a_NeCandiClu);
						reftightdis.push_back( DisMatrix_Core_Neutral[i][j][0] );
					}
					else if(DisMatrix_Core_Neutral[i][j][0] < 10 + 0.4*CluDepth || DisMatrix_Core_Neutral[i][j][1] < 20 + 0.5*CluDepth )
					{
						loosecandicluster.push_back(a_NeCandiClu);
						refloosedis.push_back( DisMatrix_Core_Neutral[i][j][0] );
					}
				}
				else if( currChargeCoreType == 0 ) // && a_recoP_ch->getEnergy() < 3 ) // && !FlagLowPtPz )	//
				{
					if(CluDepth < 20)
					{
						if(DisMatrix_Core_Neutral[i][j][1] < MinDisToNoClusterTrk)	//Tag minimal distance cluster... and see if it can be potentially linked.
						{
							MinDisToNoClusterTrk = DisMatrix_Core_Neutral[i][j][1];
						}
						if( MinDisToNoClusterTrk < 300 && abs(a_recoP_ch.getEnergy() - NeCandEn) < 1.5*a_recoP_ch.getEnergy() )	//some hard cut
						{
							tightcandicluster.push_back(a_NeCandiClu);
							reftightdis.push_back( DisMatrix_Core_Neutral[i][j][1] );
						}
					}
				}
				else
				{
					cout<<"Over balanced/Un matched/defined case: "<<a_recoP_ch.getEnergy()<<" ??? "<<currChargeCoreType<<endl; 
				}
			}
		}

		float totaltightcandiEn = 0; 
		float totalloosecandiEn = 0; 
		for(unsigned int s = 0; s < tightcandicluster.size(); s++)
		{
			totaltightcandiEn += tightcandicluster[s].getEnergy();
		}

		for(unsigned int s = 0; s < loosecandicluster.size(); s++)
		{
			totalloosecandiEn += loosecandicluster[s].getEnergy();
		}

		if( currChargeCoreType == 130 )
		{
			for(unsigned int i1 = 0; i1 < tightcandicluster.size(); i1++)
			{
				auto a_clu = tightcandicluster[i1];
				if( CurrClusterEnergy + a_clu.getEnergy() < CurrTrackEnergy + 1.0*sqrt(CurrTrackEnergy) )
				{
					mergedcluster.push_back( a_clu );
					CurrClusterEnergy += a_clu.getEnergy();
				}
			}
		}
		else if( currChargeCoreType == 131 )
		{
			for(unsigned int i1 = 0; i1 < tightcandicluster.size(); i1++)
			{
				auto  a_clu = tightcandicluster[i1];
				if( CurrClusterEnergy + a_clu.getEnergy() < CurrTrackEnergy + 1.0*sqrt(CurrTrackEnergy))
				{
					mergedcluster.push_back( a_clu );
					CurrClusterEnergy += a_clu.getEnergy();
				}
			}

			//Maybe Some ID over here?...	//layers & numbers...	//BS Case ID

			for(unsigned int i2 = 0; i2 < loosecandicluster.size(); i2++)
			{
				auto a_clu = loosecandicluster[i2];
				if( CurrClusterEnergy + a_clu.getEnergy() < CurrTrackEnergy + 1.0*sqrt(CurrTrackEnergy))       //Frags...Or some minmal hit cut
				{
					mergedcluster.push_back( a_clu );
					CurrClusterEnergy += a_clu.getEnergy();
				}
			}
		}
		else if( currChargeCoreType == 110 )
		{
			for(unsigned int i1 = 0; i1 < tightcandicluster.size(); i1++)
			{
				auto a_clu = tightcandicluster[i1];
				if(  CurrClusterEnergy + a_clu.getEnergy() < CurrTrackEnergy + 0.5*sqrt(CurrTrackEnergy))  
				{
					mergedcluster.push_back( a_clu );
					CurrClusterEnergy += a_clu.getEnergy();
				}
			}
		}
		else if( currChargeCoreType == 111 )
		{
			for(unsigned int i1 = 0; i1 < tightcandicluster.size(); i1++)
			{
				auto a_clu = tightcandicluster[i1];

				if(  CurrClusterEnergy + a_clu.getEnergy() < CurrTrackEnergy + 0.5*sqrt(CurrTrackEnergy))  
				{
					mergedcluster.push_back( a_clu );
					CurrClusterEnergy += a_clu.getEnergy();
				}
			}

			for(unsigned int i2 = 0; i2 < loosecandicluster.size(); i2++)
			{
				auto a_clu = loosecandicluster[i2];
				if( fabs(CurrClusterEnergy + a_clu.getEnergy() - CurrTrackEnergy) < fabs(CurrClusterEnergy - CurrTrackEnergy) )       //Frags...Or some minmal hit cut
				{
					mergedcluster.push_back( a_clu );
					CurrClusterEnergy += a_clu.getEnergy();
				}
			}	
		}
		else if( currChargeCoreType == 211 )	// Matched
		{
			for(unsigned int i1 = 0; i1 < tightcandicluster.size(); i1++)
			{
				auto a_clu = tightcandicluster[i1];

				if( CurrClusterEnergy + a_clu.getEnergy() < CurrTrackEnergy + 1.5*sqrt(CurrTrackEnergy))  
				{
					mergedcluster.push_back( a_clu );
					CurrClusterEnergy += a_clu.getEnergy();
				}
			}
		}
		else if( currChargeCoreType == 212)
		{
			for(unsigned int i1 = 0; i1 < tightcandicluster.size(); i1++)
			{
				auto a_clu = tightcandicluster[i1];
				if(  CurrClusterEnergy + a_clu.getEnergy() < CurrTrackEnergy + 1.5*sqrt(CurrTrackEnergy))  
				{
					mergedcluster.push_back( a_clu );
					CurrClusterEnergy += a_clu.getEnergy();
				}
			}

			for(unsigned int i2 = 0; i2 < loosecandicluster.size(); i2++)
			{
				auto a_clu = loosecandicluster[i2];
				if(  fabs(CurrClusterEnergy + a_clu.getEnergy() - CurrTrackEnergy) < fabs(CurrClusterEnergy - CurrTrackEnergy) )       //Frags...Or some minmal hit cut
				{
					mergedcluster.push_back( a_clu );
					CurrClusterEnergy += a_clu.getEnergy();
				}
			}
		}


		else if( currChargeCoreType == 0 && reftightdis.size() > 0)
		{
			float mindis = 1.0E10;
			int minindex = 0; 

			for(unsigned int i1 = 0; i1 < reftightdis.size(); i1 ++)
			{
				if(reftightdis[i1] < mindis)
				{
					mindis = reftightdis[i1];
					minindex = i1; 
				}
			}

			auto a_clu = tightcandicluster[minindex];	// Only 1? ...

			mergedcluster.push_back( a_clu );
		}
		else
		{
			cout<<"No_match, currChargeCoreType "<<currChargeCoreType<<endl; 
		}

		float CHCluEnergy = 0;

		for(int is = 0; is < int(mergedcluster.size()); is++)
		{       
			auto a_TBM_clu = mergedcluster[is];
			CHCluEnergy +=a_TBM_clu.getEnergy();
		}

		if( !( CHCluEnergy < 1 && CurrTrackEnergy > 5 ) || (MinDisToOtherTrack < 100) )	// Need to check if exist nearby larger charged cluster: maybe absorbed together and only left tiny MIP tail in the ECAL //* bool closeTonearByEnergeticChargedShower = 0  // MIP like; should also protect against energies
		{
			for(int i2 = 0; i2 < int(mergedcluster.size()); i2++)
			{
				auto a_TBM_clu = mergedcluster[i2];
				NNCTouchFlag[a_TBM_clu]	= 2; 		// can make use of this intereting flag...
			}

			float charge = a_chargedTrk.getTrackStates(0).omega/fabs(a_chargedTrk.getTrackStates(0).omega);

			auto chargeparticle = arborrecoparticleCol->create();
			chargeparticle.setCharge(charge);
			chargeparticle.addToTracks( a_chargedTrk );
			TVector3 Ptrack = Track_P3[a_chargedTrk];
			double currTrkP[3] = { Ptrack.X(), Ptrack.Y(), Ptrack.Z() };

			int flagEnergyFlow = 0;
			int currChargeCoreType2 = -99;


			auto chclustermerged = mergedclu_chCol->create();
			m_ArborToolLCIO->NaiveMergeCluConst(mergedcluster,chclustermerged);
			edm4hep::ConstCluster chclustermergedCon=chclustermerged;
			chargeparticle.addToClusters(chclustermergedCon);
			chargedclustercore_merged.push_back(chclustermerged);
			currChargeCoreType2 = m_ArborToolLCIO->ClusterFlag(chclustermerged, a_chargedTrk);

			if( currChargeCoreType2 == 130 || currChargeCoreType2 == 131 ) 
			{
				chargeparticle.setType( int(-13*charge) );
			}
			else if( currChargeCoreType2 == 110 || currChargeCoreType2 == 111 ) 
			{
				chargeparticle.setType( int(-11*charge) );
				if(CHCluEnergy > CurrTrackEnergy + 0.5*sqrt(CurrTrackEnergy) + 1)
				{
					flagEnergyFlow = 1; 
				}
			}
			else
			{
				chargeparticle.setType( Track_Type[a_chargedTrk] );
				if(CHCluEnergy > CurrTrackEnergy + 1.5*sqrt(CurrTrackEnergy) + 1)
				{
					flagEnergyFlow = 2;
				}
			}

			if( (currChargeCoreType2 != 130 && currChargeCoreType2 != 131 &&  CurrTrackEnergy > 15 && CHCluEnergy > 0.5 && CurrTrackEnergy > CHCluEnergy + 2.5*sqrt(CHCluEnergy) )  || (( currChargeCoreType2 == 130 || currChargeCoreType2 == 131 ) && CurrTrackEnergy > 100 && CHCluEnergy < 3 && chclustermerged.hits_size() < 20 ) ) 
			{
				chargeparticle.setEnergy( CHCluEnergy );
				edm4hep::Vector3f CorrMom = edm4hep::Vector3f(CHCluEnergy/CurrTrackEnergy*currTrkP[0], CHCluEnergy/CurrTrackEnergy*currTrkP[1], CHCluEnergy/CurrTrackEnergy*currTrkP[2] ); 
				chargeparticle.setMomentum( CorrMom );
			}
			else
			{
				chargeparticle.setEnergy( CurrTrackEnergy );
				chargeparticle.setMomentum(edm4hep::Vector3f( currTrkP[0],currTrkP[1],currTrkP[2]) );
			}

			if( flagEnergyFlow )
			{
				auto a_Ef_Ne_particle = arborrecoparticleCol->create();
				a_Ef_Ne_particle.setEnergy( CHCluEnergy - CurrTrackEnergy );
				TVector3 corePos = TVector3(chclustermerged.getPosition().x,chclustermerged.getPosition().y,chclustermerged.getPosition().z);
				float WFactor = (CHCluEnergy - CurrTrackEnergy)/corePos.Mag(); 
				edm4hep::Vector3f PFNEMom = edm4hep::Vector3f(WFactor*float(corePos.X()), WFactor*float(corePos.Y()), WFactor*float(corePos.Z()));
				a_Ef_Ne_particle.setMomentum(PFNEMom);
				a_Ef_Ne_particle.setMass( 0.0 );
				a_Ef_Ne_particle.setCharge( 0.0 );
				a_Ef_Ne_particle.setType(501);

				cout<<"Energy Flow Neutral Tagged "<<CHCluEnergy - CurrTrackEnergy<<endl; 
			}
			cout<<"a charged particle reconstructed with en:"<<chargeparticle.getEnergy()<<endl;

		}
		else	// push non valid tracks, etc to second iteration, as those for PreInteracting ones
		{
			SecondIterTracks.push_back(a_chargedTrk);
			cout<<"Second Iter Track Found"<<endl; 
		}	
	}

	std::vector<edm4hep::ConstCluster> BBCore; 
	BBCore.clear();
	for(int p6 = 0; p6 < NNeutralCluster; p6 ++)
	{
		auto c_clu = non_chargedclustercore[p6];
		if( NNCTouchFlag.find(c_clu) == NNCTouchFlag.end() )
		{
			BBCore.push_back(c_clu);
		} 
	}

	float NAMom[3] = {0, 0, 0};

	//Final Re-absorption
	std::vector<edm4hep::ConstCluster> NBBNeutral; 
	NBBNeutral.clear();

	for(int s = 0; s < int (BBCore.size()); s++)
	{
		auto a_clu = BBCore[s];
		TVector3 PosClu = TVector3(a_clu.getPosition().x,a_clu.getPosition().y,a_clu.getPosition().z);
		float Depth = 0; 
		Depth = DisSeedSurface(PosClu);
		float CoreEnCorr = m_ArborToolLCIO->ClusterEE(a_clu);

		if(m_ArborToolLCIO->newPhotonTag(a_clu)==1)
		{
			TVector3 BushSeedPos = TVector3(a_clu.getPosition().x,a_clu.getPosition().y,a_clu.getPosition().z);
			auto neutralparticle = arborrecoparticleCol->create();
			neutralparticle.setType(22);
			TVector3 PP = m_ArborToolLCIO->ClusterCoG(a_clu);
			NAMom[0] = CoreEnCorr*1.0/PP.Mag()*PP.X();
			NAMom[1] = CoreEnCorr*1.0/PP.Mag()*PP.Y();
			NAMom[2] = CoreEnCorr*1.0/PP.Mag()*PP.Z();
			neutralparticle.setEnergy( CoreEnCorr );
			neutralparticle.setMass( 0.0 );
			neutralparticle.setCharge( 0.0 );
			neutralparticle.setMomentum( edm4hep::Vector3f(NAMom[0],NAMom[1],NAMom[2]) );
			auto a_neclu =  mergedclu_neCol->create();
			m_ArborToolLCIO->NaiveCluConst(a_clu,a_neclu);
			a_neclu.setEnergy( CoreEnCorr );	//Reset...
			edm4hep::ConstCluster a_necluCon=a_neclu;
			neutralparticle.addToClusters(a_neclu);
		}
		else	// Distance to Charged Core > sth;
		{
			float MinDisToChCore = 1.0E9;
			float currDis = 0; 
			int NChCore = mergedclu_chCol->size();
			float closestChCluEn = 0; 			
			for(int t = 0; t < NChCore; t++)
			{
				auto a_chclu = (*mergedclu_chCol)[t];
				currDis = m_ArborToolLCIO->BushDis(a_chclu, a_clu);
				if(currDis < MinDisToChCore)
				{
					MinDisToChCore = currDis;
					closestChCluEn = a_chclu.getEnergy();	// Or the Trk En??
				}
			}
			if( MinDisToChCore > 0.4*(15 + closestChCluEn + Depth*0.01) || a_clu.getEnergy() > 2.0 )	//Joint Depth??
			{
				NBBNeutral.push_back(a_clu);
			}
		}
	}

	// Add: Neural Core Remerge & Energy Scale Recalculate, IsoHit Abso
	std::vector<edm4hep::Cluster> NBBAbs = m_ArborToolLCIO->ClusterHitAbsorbtion(NBBNeutral, IsoHits, 100); //_HitAbsCut);	// Huge??

	std::vector<float> BBAbsEn; 
	BBAbsEn.clear();

	for(unsigned s1 = 0; s1 < NBBAbs.size(); s1++)
	{
		BBAbsEn.push_back(NBBAbs[s1].getEnergy());
	}

	std::vector<int> BBAbsIndex = SortMeasure(BBAbsEn, 1);

	std::vector<edm4hep::ConstCluster > NeutronCore;
	std::vector<edm4hep::ConstCluster > NeutronFlag;
	NeutronCore.clear();
	NeutronFlag.clear();	

	for(unsigned int s2 = 0; s2 < NBBAbs.size(); s2++)	//Sort it; the first one must be a neutral core?
	{
		auto p_clu = NBBAbs[BBAbsIndex[s2]];
		float currCluEn = p_clu.getEnergy();
		std::vector<float> CluTime = m_ArborToolLCIO->ClusterTime(p_clu);
		if( (currCluEn > 1.0 || (currCluEn > 0.5 && s2 < 2) )&& CluTime[0] < 10)
		{
			NeutronCore.push_back(p_clu);
		}
		else
		{
			NeutronFlag.push_back(p_clu);
		}
	}

	std::vector<edm4hep::ConstCluster > Neutrons = m_ArborToolLCIO->ClusterAbsorbtion(NeutronCore, NeutronFlag, 200, 0.01);

	for(unsigned int s3 = 0; s3 < Neutrons.size(); s3++)
	{
		auto a_clu = Neutrons[s3];
		float CoreEnCorr = m_ArborToolLCIO->ClusterEE(a_clu);
		TVector3 PosClu = TVector3(a_clu.getPosition().x,a_clu.getPosition().y,a_clu.getPosition().z);
		float MinDisToChCore = 1.0E9;
		float RecoT0 = m_ArborToolLCIO->ClusterT0(a_clu);
		float currDis = 0; 
		int NChCore = mergedclu_chCol->size();
		for(int t = 0; t < NChCore; t++)
		{
			auto a_chclu=(*mergedclu_chCol)[t];
			currDis = m_ArborToolLCIO->BushDis(a_chclu, a_clu);
			if(currDis < MinDisToChCore)
			{
				MinDisToChCore = currDis;
			}
		}

		if( !(RecoT0>0.1 && RecoT0<1E8 && MinDisToChCore <12) )
		{
			if(m_ArborToolLCIO->newPhotonTag(a_clu)==1)
				cout<<"WARNING... Photons after neutron merge merged"<<endl; 
			auto neutralparticle = arborrecoparticleCol->create();
			neutralparticle.setType(21120);
			TVector3 PP = m_ArborToolLCIO->ClusterCoG(a_clu);

			if( RecoT0 > 0.16 && RecoT0 < 100 && MinDisToChCore > 50 && a_clu.hits_size() > 5 && abs(CoreEnCorr - 6) < 4 )
			{
				neutralparticle.setType(2112);
				float Dis = PP.Mag();
				float Beta = 1.0/(1 + 300*RecoT0/Dis);
				float PPN = 0.941/sqrt(1-Beta*Beta);
				float PPK = 0.53*PPN;
				if(PPN/CoreEnCorr < 1)
				{
					CoreEnCorr = PPN;
					neutralparticle.setType(2112);
				}
				else if(PPK/CoreEnCorr < 1 && PPN/CoreEnCorr > 1)
				{
					CoreEnCorr = 0.5*(PPN+PPK);
					neutralparticle.setType(2112);
				}
				else if(PPK/CoreEnCorr > 1)
				{
					CoreEnCorr = PPK;
					neutralparticle.setType(310);
				}
			}

			NAMom[0] = CoreEnCorr*1.0/PP.Mag()*PP.X();
			NAMom[1] = CoreEnCorr*1.0/PP.Mag()*PP.Y();
			NAMom[2] = CoreEnCorr*1.0/PP.Mag()*PP.Z();
			neutralparticle.setEnergy( CoreEnCorr );
			neutralparticle.setMass( 0.0 );
			neutralparticle.setCharge( 0.0 );
			neutralparticle.setMomentum( edm4hep::Vector3f(NAMom[0],NAMom[1],NAMom[2]) );
			auto a_neclu =  mergedclu_neCol->create();
			m_ArborToolLCIO->NaiveCluConst(a_clu,a_neclu);
			a_neclu.setEnergy( CoreEnCorr );       //Reset...
			edm4hep::ConstCluster a_necluCon=a_neclu;
			neutralparticle.addToClusters(a_necluCon);
		}
	}
	

}

StatusCode BushConnect::execute()
{
		try{
		BushConnect::Clean();	
		BushConnect::TrackSort(  );
		BushConnect::BushSelfMerge(  ); 	
		BushConnect::TagCore(  );		
		BushConnect::ParticleReco(  );
		}catch(GaudiException &e){}
	return StatusCode::SUCCESS;
}

StatusCode BushConnect::finalize()
{
	std::cout<<"Bush Connection Finished, ArborObject Formed"<<std::endl;	
	return GaudiAlgorithm::finalize();
}

