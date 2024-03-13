// /* -*- Mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */
// // Unit in code: mm, ns. 
// // NOTE: This digitialization highly matches detector geometry CRDEcalBarrel_v01.
// // TODO: read geometry info automatically.  

#include "CRDHcalDigiAlg.h" 

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

DECLARE_COMPONENT( CRDHcalDigiAlg )

CRDHcalDigiAlg::CRDHcalDigiAlg(const std::string& name, ISvcLocator* svcLoc)
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

StatusCode CRDHcalDigiAlg::initialize()
{
  if(_writeNtuple){
    std::string s_outfile = _filename;
    m_wfile = new TFile(s_outfile.c_str(), "recreate");
    t_simHit = new TTree("simHit", "simHit");
    
    t_simHit->Branch("simHit_x", &m_simHit_x);
    t_simHit->Branch("simHit_y", &m_simHit_y);
    t_simHit->Branch("simHit_z", &m_simHit_z);
    t_simHit->Branch("simHit_E", &m_simHit_E);
    t_simHit->Branch("simHit_steps", &m_simHit_steps);
    
    t_simHit->Branch("simHit_module", &m_simHit_module);
    t_simHit->Branch("simHit_stave", &m_simHit_stave);
    t_simHit->Branch("simHit_layer", &m_simHit_layer);
    t_simHit->Branch("simHit_tower", &m_simHit_tower);
    t_simHit->Branch("simHit_slice", &m_simHit_slice);
    t_simHit->Branch("simHit_cellID", &m_simHit_cellID);
  }
	std::cout<<"CRDHcalDigiAlg::m_scale="<<m_scale<<std::endl;
	m_geosvc = service<IGeomSvc>("GeomSvc");
	if ( !m_geosvc )  throw "CRDHcalDigiAlg :Failed to find GeomSvc ...";
	dd4hep::Detector* m_dd4hep = m_geosvc->lcdd();
	if ( !m_dd4hep )  throw "CRDHcalDigiAlg :Failed to get dd4hep::Detector ...";
	m_cellIDConverter = new dd4hep::rec::CellIDPositionConverter(*m_dd4hep);
   	m_decoder = m_geosvc->getDecoder(_readoutName);
	if (!m_decoder) {
		error() << "Failed to get the decoder. " << endmsg;
		return StatusCode::FAILURE;
	}

	//m_edmsvc = service<ICRDEcalSvc>("CRDEcalSvc");
	//if ( !m_edmsvc )  throw "CRDHcalDigiAlg :Failed to find CRDEcalSvc ...";

	rndm.SetSeed(_seed);
	std::cout<<"CRDHcalDigiAlg::initialize"<<std::endl;
	return GaudiAlgorithm::initialize();
}

StatusCode CRDHcalDigiAlg::execute()
{
// clock_t yyy_start, yyy_enddigi;
// yyy_start = clock(); // 记录开始时间

	if(_nEvt==0) std::cout<<"CRDHcalDigiAlg::execute Start"<<std::endl;
	std::cout<<"Processing event: "<<_nEvt<<std::endl;
   	if(_nEvt<_Nskip){ _nEvt++; return StatusCode::SUCCESS; }

	Clear();

 	const edm4hep::SimCalorimeterHitCollection* SimHitCol =  r_SimCaloCol.get();

	edm4hep::CalorimeterHitCollection* caloVec = w_DigiCaloCol.createAndPut();
	edm4hep::MCRecoCaloAssociationCollection* caloAssoVec = w_CaloAssociationCol.createAndPut();
  edm4hep::MCRecoCaloParticleAssociationCollection* caloMCPAssoVec = w_MCPCaloAssociationCol.createAndPut(); 

	if(SimHitCol == 0) 
	{
		std::cout<<"not found SimCalorimeterHitCollection"<< std::endl;
		return StatusCode::SUCCESS;
	}
  if(_Debug>=1) std::cout<<"digi, input sim hit size="<< SimHitCol->size() <<std::endl;


  for(int isim=0; isim<SimHitCol->size(); isim++){

    auto simhit = SimHitCol->at(isim);
    if(!simhit.isAvailable()) continue;
    if(simhit.getEnergy()==0) continue;

    unsigned long long id = simhit.getCellID();
    double Ehit = simhit.getEnergy();
    //Energy threshold
    if(Ehit<_MIPCali*_Eth_Mip) continue;

    //Global calibration. 
    //TODO: add more digitization terms here. 
    double Ehit_cali = Ehit*r_cali;

    //Loop contributions to get hit time and MCParticle. 
    double Thit_ave = 0.;
    double Ehit_raw = 0.;
    MCParticleToEnergyWeightMap MCPEnMap; MCPEnMap.clear();
    for(int iConb=0; iConb<simhit.contributions_size(); ++iConb){
      auto conb = simhit.getContributions(iConb);
      if(!conb.isAvailable()) continue;
      if(conb.getEnergy()==0) continue;

      Thit_ave += conb.getTime();
      
      auto mcp = conb.getParticle();
      MCPEnMap[mcp] += conb.getEnergy();
      Ehit_raw += conb.getEnergy();
    }
    Thit_ave = Thit_ave/simhit.contributions_size();
    //Create DigiHit
    auto digiHit = caloVec->create();
    digiHit.setCellID(id);
    digiHit.setEnergy(Ehit_cali);
    digiHit.setTime(Thit_ave);
    digiHit.setPosition(simhit.getPosition());

    //Create SimHit-DigiHit association. 
    auto rel = caloAssoVec->create();
    rel.setRec(digiHit);
    rel.setSim(simhit);
    rel.setWeight(1.);

    //Create DigiHit-MCParticle association.
    for(auto iter : MCPEnMap){
      auto rel_MC = caloMCPAssoVec->create();
      rel_MC.setRec(digiHit);
      rel_MC.setSim(iter.first);
      rel_MC.setWeight(iter.second/Ehit_raw);
    }

    if(_writeNtuple){
      m_simHit_x.push_back(digiHit.getPosition().x);
      m_simHit_y.push_back(digiHit.getPosition().y);
      m_simHit_z.push_back(digiHit.getPosition().z);
      m_simHit_E.push_back(digiHit.getEnergy());
      m_simHit_steps.push_back(simhit.contributions_size());
      m_simHit_module.push_back(m_decoder->get(id, "module"));
      m_simHit_stave.push_back(m_decoder->get(id, "stave"));
      m_simHit_layer.push_back(m_decoder->get(id, "layer"));
      m_simHit_slice.push_back(m_decoder->get(id, "slice"));
      m_simHit_tower.push_back(m_decoder->get(id, "tower"));
      m_simHit_cellID.push_back(id);
    }
  }

	if(_writeNtuple) t_simHit->Fill();

	_nEvt ++ ;
	return StatusCode::SUCCESS;
}

StatusCode CRDHcalDigiAlg::finalize()
{
  if(_writeNtuple){
	  m_wfile->cd();
	  t_simHit->Write();
    m_wfile->Close();
	  delete m_wfile, t_simHit; 
  }

	info() << "Processed " << _nEvt << " events " << endmsg;
	delete m_cellIDConverter, m_decoder, m_geosvc;
	return GaudiAlgorithm::finalize();
}


void CRDHcalDigiAlg::Clear(){
	m_simHit_x.clear();
	m_simHit_y.clear();
	m_simHit_z.clear();
	m_simHit_E.clear();
	m_simHit_steps.clear();
	m_simHit_module.clear();
	m_simHit_stave.clear();
	m_simHit_layer.clear();
	m_simHit_slice.clear();
	m_simHit_tower.clear();
  m_simHit_cellID.clear();
}

