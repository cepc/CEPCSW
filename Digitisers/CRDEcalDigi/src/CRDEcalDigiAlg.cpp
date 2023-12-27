/* -*- Mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */
// Unit in code: mm, ns. 
// NOTE: This digitialization highly matches detector geometry CRDEcalBarrel_v01.
// TODO: read geometry info automatically.  
#include "CRDEcalDigiAlg.h"

#include "edm4hep/SimCalorimeterHit.h"
#include "edm4hep/CalorimeterHit.h"
#include "edm4hep/Vector3f.h"
#include "edm4hep/Cluster.h"

#include "DD4hep/Detector.h"
#include <DD4hep/Objects.h>
#include <DDRec/CellIDPositionConverter.h>

#include "TVector3.h"
#include <math.h>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <map>

// #include <fstream>
// #include <ctime>

#define C 299.79  // unit: mm/ns
#define PI 3.141592653
using namespace std;
using namespace dd4hep;

DECLARE_COMPONENT( CRDEcalDigiAlg )

CRDEcalDigiAlg::CRDEcalDigiAlg(const std::string& name, ISvcLocator* svcLoc)
  : GaudiAlgorithm(name, svcLoc),
    _nEvt(0)
{
  
	// Input collections
	declareProperty("SimCaloHitCollection", r_SimCaloCol, "Handle of the Input SimCaloHit collection");
  
	// Output collections
	declareProperty("CaloHitCollection", w_DigiCaloCol, "Handle of Digi CaloHit collection");
	declareProperty("CaloAssociationCollection", w_CaloAssociationCol, "Handle of CaloAssociation collection");
	declareProperty("CaloMCPAssociationCollection", w_MCPCaloAssociationCol, "Handle of CaloAssociation collection");
   
}

StatusCode CRDEcalDigiAlg::initialize()
{

  if(_writeNtuple){
    std::string s_outfile = _filename;
    m_wfile = new TFile(s_outfile.c_str(), "recreate");
    t_SimCont = new TTree("SimStep", "SimStep");
    t_SimBar = new TTree("SimBarHit", "SimBarHit");
    t_SimCont->Branch("step_x", &m_step_x);
    t_SimCont->Branch("step_y", &m_step_y);
    t_SimCont->Branch("step_z", &m_step_z);
    t_SimCont->Branch("step_t", &m_step_t);			// yyy: time of each step
    t_SimCont->Branch("stepBar_x", &m_stepBar_x);
    t_SimCont->Branch("stepBar_y", &m_stepBar_y);
    t_SimCont->Branch("stepBar_z", &m_stepBar_z);
    t_SimCont->Branch("step_E", &m_step_E);
    t_SimCont->Branch("step_T1", &m_step_T1);
    t_SimCont->Branch("step_T2", &m_step_T2);
    t_SimBar->Branch("totE", &totE);
    t_SimBar->Branch("simBar_x", &m_simBar_x);
    t_SimBar->Branch("simBar_y", &m_simBar_y);
    t_SimBar->Branch("simBar_z", &m_simBar_z);
    t_SimBar->Branch("simBar_T1", &m_simBar_T1);
    t_SimBar->Branch("simBar_T2", &m_simBar_T2);
    t_SimBar->Branch("simBar_Q1", &m_simBar_Q1);
    t_SimBar->Branch("simBar_Q2", &m_simBar_Q2);
    t_SimBar->Branch("simBar_module", &m_simBar_module);
    t_SimBar->Branch("simBar_stave", &m_simBar_stave);
    t_SimBar->Branch("simBar_dlayer", &m_simBar_dlayer);
    t_SimBar->Branch("simBar_part", &m_simBar_part);
    t_SimBar->Branch("simBar_slayer", &m_simBar_slayer);
    t_SimBar->Branch("simBar_cellID", &m_simBar_cellID);
  }

	std::cout<<"CRDEcalDigiAlg::m_scale="<<m_scale<<std::endl;
	m_geosvc = service<IGeomSvc>("GeomSvc");
	if ( !m_geosvc )  throw "CRDEcalDigiAlg :Failed to find GeomSvc ...";
	dd4hep::Detector* m_dd4hep = m_geosvc->lcdd();
	if ( !m_dd4hep )  throw "CRDEcalDigiAlg :Failed to get dd4hep::Detector ...";
	m_cellIDConverter = new dd4hep::rec::CellIDPositionConverter(*m_dd4hep);
   m_decoder = m_geosvc->getDecoder(_readoutName);
	if (!m_decoder) {
		error() << "Failed to get the decoder. " << endmsg;
		return StatusCode::FAILURE;
	}

	//m_edmsvc = service<ICRDEcalSvc>("CRDEcalSvc");
	//if ( !m_edmsvc )  throw "CRDEcalDigiAlg :Failed to find CRDEcalSvc ...";

	rndm.SetSeed(_seed);
	std::cout<<"CRDEcalDigiAlg::initialize"<<std::endl;
	return GaudiAlgorithm::initialize();
}

StatusCode CRDEcalDigiAlg::execute()
{
// clock_t yyy_start, yyy_enddigi;
// yyy_start = clock(); // 记录开始时间

	if(_nEvt==0) std::cout<<"CRDEcalDigiAlg::execute Start"<<std::endl;
	std::cout<<"Processing event: "<<_nEvt<<std::endl;
   if(_nEvt<_Nskip){ _nEvt++; return StatusCode::SUCCESS; }

	Clear();

 	const edm4hep::SimCalorimeterHitCollection* SimHitCol =  r_SimCaloCol.get();

	edm4hep::CalorimeterHitCollection* caloVec = w_DigiCaloCol.createAndPut();
	edm4hep::MCRecoCaloAssociationCollection* caloAssoVec = w_CaloAssociationCol.createAndPut();
  edm4hep::MCRecoCaloParticleAssociationCollection* caloMCPAssoVec = w_MCPCaloAssociationCol.createAndPut();
 	std::vector<edm4hep::SimCalorimeterHit> m_simhitCol; m_simhitCol.clear();
  std::vector<CaloBar> m_barCol; m_barCol.clear(); 

  if(SimHitCol == 0) 
  {
     std::cout<<"not found SimCalorimeterHitCollection"<< std::endl;
     return StatusCode::SUCCESS;
  }
  if(_Debug>=1) std::cout<<"digi, input sim hit size="<< SimHitCol->size() <<std::endl;

	totE=0;
	double totE_Digi=0;

	//Merge input simhit(steps) to real simhit(bar).
	MergeHits(*SimHitCol, m_simhitCol);
	if(_Debug>=1) std::cout<<"Finish Hit Merge, with Nhit: "<<m_simhitCol.size()<<std::endl;

	//Loop in SimHit, digitalize SimHit to DigiBar
	for(int i=0;i<m_simhitCol.size();i++){

		auto SimHit = m_simhitCol.at(i);
		if(SimHit.getEnergy()<_Eth) continue;

		unsigned long long id = SimHit.getCellID();
		CaloBar hitbar;
		hitbar.setcellID( id);
		hitbar.setcellID(	m_decoder->get(id, "system"), 
											m_decoder->get(id, "module"), 
											m_decoder->get(id, "stave"), 
											m_decoder->get(id, "dlayer"), 
											m_decoder->get(id, "part"), 
											m_decoder->get(id, "slayer"),
											m_decoder->get(id, "bar"));

		double Lbar = GetBarLength(hitbar);  //NOTE: Is fixed with geometry CRDEcalBarrel_v01. 
		dd4hep::Position hitpos = m_cellIDConverter->position(id);
    TVector3 barpos(10*hitpos.x(), 10*hitpos.y(), 10*hitpos.z()); //cm to mm.
		hitbar.setPosition(barpos);
		//hitbar.T1 = 99999; hitbar.T2 = 99999;
		//if(_Debug>=2) std::cout<<"SimHit contribution size: "<<SimHit.contributions_size()<<std::endl;

    MCParticleToEnergyWeightMap MCPEnMap; MCPEnMap.clear();
		std::vector<HitStep> DigiLvec; DigiLvec.clear();
		std::vector<HitStep> DigiRvec; DigiRvec.clear();
		double totQ1 = 0;
		double totQ2 = 0;

		//Loop in all SimHitContribution(G4Step). 
		for(int iCont=0; iCont < SimHit.contributions_size(); ++iCont){
			auto conb = SimHit.getContributions(iCont);
			if( !conb.isAvailable() ) { std::cout<<"CRDEcalDigiAlg  Can not get SimHitContribution: "<<iCont<<std::endl; continue;}

			double en = conb.getEnergy();
			if(en == 0) continue;

      auto mcp = conb.getParticle();
      MCPEnMap[mcp] += en;
			TVector3 steppos(conb.getStepPosition().x, conb.getStepPosition().y, conb.getStepPosition().z);
			TVector3 rpos = steppos-hitbar.getPosition();
			float step_time = conb.getTime();		// yyy: step time

			m_step_x.push_back(steppos.x());
			m_step_y.push_back(steppos.y());
			m_step_z.push_back(steppos.z());
			m_step_t.push_back(step_time);			// yyy: push back step time
			m_step_E.push_back(en);
			m_stepBar_x.push_back(hitbar.getPosition().x());
			m_stepBar_y.push_back(hitbar.getPosition().y());
			m_stepBar_z.push_back(hitbar.getPosition().z());

			if(_Debug>=3){
				cout<<"Cell Pos: "<<hitbar.getPosition().x()<<'\t'<<hitbar.getPosition().y()<<'\t'<<hitbar.getPosition().z()<<endl;
				cout<<"step pos: "<<steppos.x()<<'\t'<<steppos.y()<<'\t'<<steppos.z()<<endl;
				cout<<"Relative pos: "<<rpos.x()<<'\t'<<rpos.y()<<'\t'<<rpos.z()<<endl;
				cout<<"Cell: "<<hitbar.getModule()<<"  "<<hitbar.getDlayer()<<"  "<<hitbar.getSlayer()<<endl;
			}

			//Get digitalized signal(Q1, Q2, T1, T2) from step
			//Define: 1 is left, 2 is right, clockwise direction in phi. 

			int sign=-999;
			if(hitbar.getSlayer()==1) sign = rpos.z()==0 ? 1 : rpos.z()/fabs(rpos.z());
			else{
				if(hitbar.getModule()==0 || hitbar.getModule()==1 || hitbar.getModule()==7) sign = rpos.x()==0 ?  1: rpos.x()/fabs(rpos.x());
				if(hitbar.getModule()==3 || hitbar.getModule()==4 || hitbar.getModule()==5) sign = rpos.x()==0 ? -1:-rpos.x()/fabs(rpos.x());
				else if(hitbar.getModule()==2) sign = rpos.y()==0 ?  1: rpos.y()/fabs(rpos.y());
				else if(hitbar.getModule()==6) sign = rpos.y()==0 ? -1:-rpos.y()/fabs(rpos.y());
			}
			if(!fabs(sign)) {std::cout<<"ERROR: Wrong bar direction/position!"<<std::endl; continue;}


			double Qi_left = en*exp(-(Lbar/2 + sign*sqrt(rpos.Mag2()))/Latt);	
			double Qi_right = en*exp(-(Lbar/2 - sign*sqrt(rpos.Mag2()))/Latt);

			if(_Debug>=3){
				cout<<Qi_left<<'\t'<<Qi_right<<endl;
				cout<<Lbar<<'\t'<<sign*sqrt(rpos.Mag2())<<endl;
			}


			double Ti_left = -1; int looptime=0;
			while(Ti_left<0){ 
				// Ti_left = Tinit + rndm.Gaus(nMat*(Lbar/2 + sign*sqrt(rpos.Mag2()))/C, Tres); 
				Ti_left = Tinit + rndm.Gaus(nMat*(Lbar/2 + sign*sqrt(rpos.Mag2()))/C, Tres) + step_time;  // yyy: add step time 
				looptime++;
				if(looptime>500){ std::cout<<"ERROR: Step "<<iCont<<" can not get a positive left-side time!"<<std::endl; break;}
			}
			if(looptime>500) continue;		
			double Ti_right = -1; looptime=0;
			while(Ti_right<0){ 
				// Ti_right = Tinit + rndm.Gaus(nMat*(Lbar/2 - sign*sqrt(rpos.Mag2()))/C, Tres); 
				Ti_right = Tinit + rndm.Gaus(nMat*(Lbar/2 - sign*sqrt(rpos.Mag2()))/C, Tres) + step_time;  // yyy: add step time 
				looptime++;
            if(looptime>500){ std::cout<<"ERROR: Step "<<iCont<<" can not get a positive right-side time!"<<std::endl; break;}
			}
			if(looptime>500) continue;		


			m_step_T1.push_back(Ti_left);
			m_step_T2.push_back(Ti_right);
			totQ1 += Qi_left;
			totQ2 += Qi_right;
	
			HitStep stepoutL, stepoutR;
			stepoutL.setQ(Qi_left); stepoutL.setT(Ti_left);
			stepoutR.setQ(Qi_right); stepoutR.setT(Ti_right);
			DigiLvec.push_back(stepoutL);
			DigiRvec.push_back(stepoutR);

		}

		//Time digitalization
		//if(_Debug>=2) std::cout<<"Time Digitalize: time at Q >"<<_Qthfrac<<"*totQ"<<std::endl;
		std::sort(DigiLvec.begin(), DigiLvec.end());
		std::sort(DigiRvec.begin(), DigiRvec.end());
		double thQ1=0;
		double thQ2=0;
		double thT1, thT2; 
		for(int iCont=0;iCont<DigiLvec.size();iCont++){
			thQ1 += DigiLvec[iCont].getQ();
			if(thQ1>totQ1*_Qthfrac){ 
				thT1 = DigiLvec[iCont].getT(); 
				if(_Debug>=3) std::cout<<"Get T1 at index: "<<iCont<<std::endl;
				break;
			}
		}
		for(int iCont=0;iCont<DigiRvec.size();iCont++){
			thQ2 += DigiRvec[iCont].getQ();
			if(thQ2>totQ2*_Qthfrac){ 
				thT2 = DigiRvec[iCont].getT(); 
				if(_Debug>=3) std::cout<<"Get T2 at index: "<<iCont<<std::endl;
				break;
			}
		}
		hitbar.setQ(totQ1, totQ2);
		hitbar.setT(thT1, thT2);
    //End bar digitization. 

    //2 hits with double-readout time. 
    edm4hep::Vector3f m_pos(hitbar.getPosition().X(), hitbar.getPosition().Y(), hitbar.getPosition().Z());
    auto digiHit1 = caloVec->create();
    digiHit1.setCellID(hitbar.getcellID());
    digiHit1.setEnergy(hitbar.getQ1());
    digiHit1.setTime(hitbar.getT1());
    digiHit1.setPosition(m_pos);
    auto digiHit2 = caloVec->create();
    digiHit2.setCellID(hitbar.getcellID());
    digiHit2.setEnergy(hitbar.getQ2());
    digiHit2.setTime(hitbar.getT2());
    digiHit2.setPosition(m_pos);

    //SimHit - CaloHit association
    auto rel1 = caloAssoVec->create();
    rel1.setRec(digiHit1);
    rel1.setSim(SimHit);
    rel1.setWeight( hitbar.getQ1()/(hitbar.getQ1()+hitbar.getQ2()) );
    auto rel2 = caloAssoVec->create();
    rel2.setRec(digiHit2);
    rel2.setSim(SimHit);
    rel2.setWeight( hitbar.getQ2()/(hitbar.getQ1()+hitbar.getQ2()) );


    //MCParticle - CaloHit association
    //float maxMCE = -99.;
    //edm4hep::MCParticle selMCP; 
    for(auto iter : MCPEnMap){
      //if(iter.second>maxMCE){
      //  maxMCE = iter.second;
      //  selMCP = iter.first;
      //}
      auto rel_MCP1 = caloMCPAssoVec->create();
      rel_MCP1.setRec(digiHit1);
      rel_MCP1.setSim(iter.first);
      rel_MCP1.setWeight(iter.second/digiHit1.getEnergy());
      auto rel_MCP2 = caloMCPAssoVec->create();
      rel_MCP2.setRec(digiHit2);
      rel_MCP2.setSim(iter.first);
      rel_MCP2.setWeight(iter.second/digiHit2.getEnergy());      
    }

    //if(selMCP.isAvailable()){
    //  auto rel_MCP1 = caloMCPAssoVec->create();
    //  rel_MCP1.setRec(digiHit1);
    //  rel_MCP1.setSim(selMCP);
    //  rel_MCP1.setWeight(1.);
    //  auto rel_MCP2 = caloMCPAssoVec->create();
    //  rel_MCP2.setRec(digiHit2);
    //  rel_MCP2.setSim(selMCP);
    //  rel_MCP2.setWeight(1.);
    //}

    m_barCol.push_back(hitbar);
		totE+=(hitbar.getQ1()+hitbar.getQ2())/2;
		
    //Temp: write into trees. 
    if(_writeNtuple){
      m_simBar_x.push_back(hitbar.getPosition().x());
      m_simBar_y.push_back(hitbar.getPosition().y());
      m_simBar_z.push_back(hitbar.getPosition().z());
      m_simBar_Q1.push_back(hitbar.getQ1());
      m_simBar_Q2.push_back(hitbar.getQ2());
      m_simBar_T1.push_back(hitbar.getT1());
      m_simBar_T2.push_back(hitbar.getT2());
      m_simBar_module.push_back(hitbar.getModule());
      m_simBar_stave.push_back(hitbar.getStave());
      m_simBar_dlayer.push_back(hitbar.getDlayer());
      m_simBar_part.push_back(hitbar.getPart());
      m_simBar_slayer.push_back(hitbar.getSlayer());
      m_simBar_cellID.push_back(hitbar.getcellID());
    }
	}


  if(_writeNtuple){
  	t_SimCont->Fill();
	  t_SimBar->Fill();
  }
	if(_Debug>=1) std::cout<<"End Loop: Bar Digitalization!"<<std::endl;
	std::cout<<"Total Bar Energy: "<<totE<<std::endl;

// yyy_enddigi = clock();
// double duration_digi = double(yyy_enddigi - yyy_start) / CLOCKS_PER_SEC;
// // 将时间输出到txt文件中
// std::ofstream outfile("runtime_ecaldigi.txt", std::ios::app);
// outfile << _nEvt << "    " << duration_digi << std::endl;
// outfile.close();

  _nEvt ++ ;
  //delete SimHitCol, caloVec, caloAssoVec; 
  m_simhitCol.clear();
  return StatusCode::SUCCESS;
}

StatusCode CRDEcalDigiAlg::finalize()
{
  if(_writeNtuple){
  	m_wfile->cd();
	  t_SimCont->Write();
  	t_SimBar->Write();
	  m_wfile->Close();
    delete m_wfile, t_SimCont, t_SimBar; 
  }

  info() << "Processed " << _nEvt << " events " << endmsg;
  delete m_cellIDConverter, m_decoder, m_geosvc;
  return GaudiAlgorithm::finalize();
}

StatusCode CRDEcalDigiAlg::MergeHits( const edm4hep::SimCalorimeterHitCollection& m_col, std::vector<edm4hep::SimCalorimeterHit>& m_hits ){

  m_hits.clear(); 
	std::vector<edm4hep::MutableSimCalorimeterHit> m_mergedhit;
	m_mergedhit.clear();

	for(int iter=0; iter<m_col.size(); iter++){
		edm4hep::SimCalorimeterHit m_step = m_col[iter];
		if(!m_step.isAvailable()){ cout<<"ERROR HIT!"<<endl; continue;}
		if(m_step.getEnergy()==0 || m_step.contributions_size()<1) continue;
		unsigned long long cellid = m_step.getCellID();
		dd4hep::Position hitpos = m_cellIDConverter->position(cellid);
		edm4hep::Vector3f pos(hitpos.x()*10, hitpos.y()*10, hitpos.z()*10);

		edm4hep::MutableCaloHitContribution conb;
		conb.setEnergy(m_step.getEnergy());
		conb.setStepPosition(m_step.getPosition());
    conb.setParticle( m_step.getContributions(0).getParticle() );
		conb.setTime(m_step.getContributions(0).getTime());

		edm4hep::MutableSimCalorimeterHit m_hit = find(m_mergedhit, cellid);
		if(m_hit.getCellID()==0){
			//m_hit = new edm4hep::SimCalorimeterHit();
			m_hit.setCellID(cellid);
			m_hit.setPosition(pos);
			m_mergedhit.push_back(m_hit);
		}
		m_hit.addToContributions(conb);
		m_hit.setEnergy(m_hit.getEnergy()+m_step.getEnergy());
	}

  for(auto iter = m_mergedhit.begin(); iter!=m_mergedhit.end(); iter++){
    edm4hep::SimCalorimeterHit constsimhit = *iter; 
    m_hits.push_back( constsimhit );  
  }
  return StatusCode::SUCCESS; 
}



double CRDEcalDigiAlg::GetBarLength(CaloBar& bar){
	//TODO: reading bar length from geosvc. 
	if(bar.getSlayer()==1) return 600.;
	else return 480.-bar.getDlayer()*10.;
}

/*
dd4hep::Position CRDEcalDigiAlg::GetCellPos(dd4hep::Position& pos, CaloBar& bar){
	dd4hep::Position rpos = pos-bar.getPosition();
	TVector3 vec(0,0,0); 
	if(bar.getSlayer()==1) vec.SetXYZ(0, 0, floor(rpos.z()/10)*10+5 );
	else if(bar.getSlayer()==0){
		if((bar.getModule()==0||bar.getModule()==4) && bar.getDlayer()%2==1) vec.SetXYZ(floor(rpos.x()/10)*10+5,0,0);
		if((bar.getModule()==0||bar.getModule()==4) && bar.getDlayer()%2==0) vec.SetXYZ(floor((rpos.x()-5)/10)*10+10,0,0);
		if((bar.getModule()==2||bar.getModule()==6) && bar.getDlayer()%2==1) vec.SetXYZ(0, floor(rpos.y()/10)*10+5,0);
		if((bar.getModule()==2||bar.getModule()==6) && bar.getDlayer()%2==0) vec.SetXYZ(0, floor((rpos.y()-5)/10)*10+10,0);
		if(bar.getModule()==1 || bar.getModule()==5){
			TVector3 unitv(1./sqrt(2), -1./sqrt(2), 0);
			if(bar.getDlayer()%2==1) vec = (floor(rpos.Dot(unitv)/10)*10+5)*unitv;
			if(bar.getDlayer()%2==0) vec = (floor((rpos.Dot(unitv)-5)/10)*10+10)*unitv;
		}
		if(bar.getModule()==3 || bar.getModule()==7){
			TVector3 unitv(1./sqrt(2), 1./sqrt(2), 0);
			if(bar.getDlayer()%2==1) vec = (floor(rpos.Dot(unitv)/10)*10+5)*unitv;
			if(bar.getDlayer()%2==0) vec = (floor((rpos.Dot(unitv)-5)/10)*10+10)*unitv;
		}
	}
	dd4hep::Position relv(vec.x(), vec.y(), vec.z());
	return relv+bar.getPosition();
}


edm4hep::MutableSimCalorimeterHit CRDEcalDigiAlg::find(edm4hep::SimCalorimeterHitCollection& m_col, dd4hep::Position& pos){
   for(int i=0;i<m_col.size();i++){
    edm4hep::MutableSimCalorimeterHit hit = m_col[i];
		dd4hep::Position ipos(hit.getPosition().x, hit.getPosition().y, hit.getPosition().z);
		if(ipos==pos) return hit;
	}
   edm4hep::MutableSimCalorimeterHit hit;
   hit.setCellID(0);
   return hit;
}
*/
edm4hep::MutableSimCalorimeterHit CRDEcalDigiAlg::find(const std::vector<edm4hep::MutableSimCalorimeterHit>& m_col, unsigned long long& cellid) const{
   for(int i=0;i<m_col.size();i++){
		edm4hep::MutableSimCalorimeterHit hit=m_col.at(i);
		if(hit.getCellID() == cellid) return hit;
	}
	edm4hep::MutableSimCalorimeterHit hit ;
	hit.setCellID(0);
	return hit;
}

void CRDEcalDigiAlg::Clear(){
  totE = -99;
	m_step_x.clear();
	m_step_y.clear();
	m_step_z.clear();
	m_step_t.clear();   // yyy: clear
	m_step_E.clear();
	m_stepBar_x.clear();
	m_stepBar_y.clear();
	m_stepBar_z.clear();
	m_step_T1.clear();
	m_step_T2.clear();
	m_simBar_x.clear();
	m_simBar_y.clear();
	m_simBar_z.clear();
	m_simBar_T1.clear();
	m_simBar_T2.clear();
	m_simBar_Q1.clear();
	m_simBar_Q2.clear();
	m_simBar_module.clear();
	m_simBar_stave.clear();
	m_simBar_dlayer.clear();
	m_simBar_part.clear();
	m_simBar_slayer.clear();
  m_simBar_cellID.clear();
}

