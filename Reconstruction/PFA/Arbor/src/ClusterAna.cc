#include "ClusterAna.hh"
#include <EVENT/LCCollection.h>
#include <IMPL/LCCollectionVec.h>
#include <EVENT/LCFloatVec.h>
#include <EVENT/MCParticle.h>
#include <EVENT/ReconstructedParticle.h>
#include <IMPL/MCParticleImpl.h>
#include <values.h>
#include <string>
#include <iostream>
#include <EVENT/LCFloatVec.h>
#include <EVENT/LCParameters.h>
#include <stdexcept>
#include <TFile.h> 
#include <TTree.h>
#include <TH1F.h>
#include <TVector3.h>
#include <TRandom.h>
#include <Rtypes.h> 
#include <sstream>		
#include <cmath>
#include <vector>
#include <TMath.h>
#include "TLorentzVector.h"
#include <UTIL/CellIDDecoder.h>

using namespace std;

DECLARE_COMPONENT(ClusterAna)


const float sqrts = 240.0; 	//GeV

const string ECALCellIDDecoder = "M:3,S-1:3,I:9,J:9,K-1:6";
//TH1F *h_hit;

ClusterAna::ClusterAna(const std::string& name, ISvcLocator* svcLoc)
    : GaudiAlgorithm(name, svcLoc),
      _output(0)
{
}

StatusCode ClusterAna::initialize() {

    // printParameters();

    TFile *tree_file=new TFile(_treeFileName.value().c_str(),(_overwrite ? "RECREATE" : "UPDATE"));

    if (!tree_file->IsOpen()) {
        delete tree_file;
        tree_file=new TFile(_treeFileName.value().c_str(),"NEW");
    }

    //h_hit=new TH1F("hit","hit",80,0,80);
    _outputTree = new TTree(_treeName.value().c_str(),_treeName.value().c_str());
    _outputTree->SetAutoSave(32*1024*1024);  // autosave every 32MB
    _outputTree->Branch("EventNr", &_eventNr, "EventNr/I");
    _outputTree->Branch("Num", &_Num, "Num/I");
    _outputTree->Branch("NClu", &_NClu, "NClu/I");
    _outputTree->Branch("NPFO", &_NPFO, "NPFO/I");
    _outputTree->Branch("EcalTotalE", &_EcalTotalE, "EcalTotalE/F");
    _outputTree->Branch("HcalTotalE", &_HcalTotalE, "HcalTotalE/F");
    _outputTree->Branch("CluE", &_CluE, "CluE/F");

	_outputClu = new TTree("Clu","Clu");
    _outputClu->Branch("EventNr", &_eventNr, "EventNr/I");
    _outputClu->Branch("Num", &_Num, "Num/I");
    _outputClu->Branch("Nhit", &_Nhit, "Nhit/I");
    _outputClu->Branch("CluEn", &_CluEn, "CluEn/F");
    _outputClu->Branch("CluPos", CluPos, "CluPos[3]/F");

    _outputPFO = new TTree("PFO","PFO");
    _outputPFO->Branch("Num", &_Num, "Num/I");
    _outputPFO->Branch("PID", &_PID, "PID/I");
    _outputPFO->Branch("PFOEn", &_PFOEn, "PFOEn/F");


    _Num = 0;

     return GaudiAlgorithm::initialize();

}

StatusCode ClusterAna::execute() 
{		

    EVENT::LCEvent* evtP = nullptr;

    std::vector<CaloHitColHandler*> hdl_EcalHitColl{
        &m_ecalbarrelhitcol,
            &m_ecalendcaphitcol
            };
    std::vector<CaloHitColHandler*> hdl_HcalHitColl{
        &m_hcalbarrelhitcol,
            &m_hcalendcaphitcol,
            &m_hcalotherhitcol
            };

    std::vector<std::string> EcalHitColl;
    std::vector<std::string> HcalHitColl;
    EcalHitColl.push_back("ECALBarrel");
    EcalHitColl.push_back("ECALEndcap");
    HcalHitColl.push_back("HCALBarrel");
    HcalHitColl.push_back("HCALEndcap");
    HcalHitColl.push_back("HCALOther");
    _CluE=0;
	
    try{

        for(int t = 0; t< int(hdl_EcalHitColl.size()); t++) {
            const edm4hep::CalorimeterHitCollection* ecalcoll = hdl_EcalHitColl[t]->get();
            for(auto hit: *ecalcoll) {
                // TODO
                int NLayer = 0;
                _EcalTotalE += hit.getEnergy();

            }
        }

//        for(int t2 = 0; t2< int(hdl_HcalHitColl.size()); t2++) {
//            const edm4hep::CalorimeterHitCollection* hcalcoll = hdl_HcalHitColl[t2]->get();
//            for (auto hit: *hcalcoll) {
//                // TODO
//                int NLayer = 0;
//                int HLayer = NLayer+30;
//                // UTIL::CellIDDecoder<EVENT::CalorimeterHit> idDecoder(ECALCellIDDecoder);
//                // int NLayer = idDecoder(a_hit)["K-1"];
//
//                //h_hit->Fill(HLayer,a_hit->getEnergy());
//                _HcalTotalE += hit.getEnergy();
//                
//            }
//        }
//
//
    }catch(lcio::DataNotAvailableException err) { }

    try{
        auto col_Clu = m_Clu.get();
	_NClu=col_Clu->size();
        for(int i0 = 0; i0 < col_Clu->size(); i0++) {
            auto a_clu = (*col_Clu)[i0];
            auto currPos0 = a_clu.getPosition();
            TVector3 currPos(currPos0.x, currPos0.y, currPos0.z);
	    _CluEn=a_clu.getEnergy();
	    _CluE+=a_clu.getEnergy();
            

            _outputClu->Fill();

        }
    }catch (lcio::DataNotAvailableException err) { }

    try{
	    auto col_PFO =  m_PFO.get();
	    _NPFO=col_PFO->size();
	    for(int i0 = 0; i0 < col_PFO->size(); i0++){
		    auto a_PFO=(*col_PFO)[i0];
		    _PFOEn=a_PFO.getEnergy();
		    _PID=a_PFO.getType();
		    _outputPFO->Fill();
	    }
    }catch (lcio::DataNotAvailableException err) { }
    _outputTree->Fill();
    _Num++;
    // }  	  
    return StatusCode::SUCCESS;

}	

StatusCode ClusterAna::finalize()
{

    if (_outputTree) {

        TFile *tree_file = _outputTree->GetCurrentFile(); //just in case we switched to a new file
        tree_file->Write();
        delete tree_file;
    }

    return GaudiAlgorithm::finalize();
}



