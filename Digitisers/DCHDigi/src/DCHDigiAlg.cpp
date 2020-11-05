/* -*- Mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */
#include "DCHDigiAlg.h"


#include "edm4hep/SimCalorimeterHit.h"
#include "edm4hep/CalorimeterHit.h"
#include "edm4hep/Vector3f.h"

#include "DD4hep/Detector.h"
#include <DD4hep/Objects.h>
#include "DDRec/Vector3D.h"


#include <array>
#include <math.h>
#include <cmath>
#include <algorithm>

DECLARE_COMPONENT( DCHDigiAlg )

DCHDigiAlg::DCHDigiAlg(const std::string& name, ISvcLocator* svcLoc)
  : GaudiAlgorithm(name, svcLoc),
    _nEvt(0)
{
  
  // Input collections
  declareProperty("SimDCHitCollection", r_SimDCHCol, "Handle of the Input SimHit collection");
  
  // Output collections
  declareProperty("DigiDCHitCollection", w_DigiDCHCol, "Handle of Digi DCHit collection");
  
  declareProperty("AssociationCollection", w_AssociationCol, "Handle of Association collection");
   
}

StatusCode DCHDigiAlg::initialize()
{
  m_geosvc = service<IGeomSvc>("GeomSvc");
  if ( !m_geosvc )  throw "DCHDigiAlg :Failed to find GeomSvc ...";
  dd4hep::Detector* m_dd4hep = m_geosvc->lcdd();
  if ( !m_dd4hep )  throw "DCHDigiAlg :Failed to get dd4hep::Detector ...";
  dd4hep::Readout readout = m_dd4hep->readout(m_readout_name);
  m_segmentation = dynamic_cast<dd4hep::DDSegmentation::GridDriftChamber*>(readout.segmentation().segmentation());
  m_decoder = m_geosvc->getDecoder(m_readout_name);
  if (!m_decoder) {
      error() << "Failed to get the decoder. " << endmsg;
      return StatusCode::FAILURE;
  }

  std::string s_output=m_Output; 
  if(m_WriteAna){
      m_fout = new TFile(s_output.c_str(),"RECREATE"); 
      m_tree = new TTree("evt","tree");
      m_tree->Branch("chamber"  , &m_chamber  );
      m_tree->Branch("layer"    , &m_layer    );
      m_tree->Branch("cell"     , &m_cell     );
      m_tree->Branch("cell_x"   , &m_cell_x   );
      m_tree->Branch("cell_y"   , &m_cell_y   );
      m_tree->Branch("simhit_x" , &m_simhit_x );
      m_tree->Branch("simhit_y" , &m_simhit_y );
      m_tree->Branch("simhit_z" , &m_simhit_z );
      m_tree->Branch("hit_x"    , &m_hit_x    );
      m_tree->Branch("hit_y"    , &m_hit_y    );
      m_tree->Branch("hit_z"    , &m_hit_z    );
      m_tree->Branch("dca"      , &m_dca      );
      m_tree->Branch("hit_dE"   , &m_hit_dE   );
      m_tree->Branch("hit_dE_dx", &m_hit_dE_dx);
  }
  std::cout<<"DCHDigiAlg::initialized"<< std::endl;
  return GaudiAlgorithm::initialize();
}

void DCHDigiAlg::Reset()
{

  std::vector<int  >().swap( m_chamber   );
  std::vector<int  >().swap( m_layer     );
  std::vector<int  >().swap( m_cell      );
  std::vector<float>().swap( m_cell_x    );
  std::vector<float>().swap( m_cell_y    );
  std::vector<float>().swap( m_simhit_x  );
  std::vector<float>().swap( m_simhit_y  );
  std::vector<float>().swap( m_simhit_z  );
  std::vector<float>().swap( m_hit_x     );
  std::vector<float>().swap( m_hit_y     );
  std::vector<float>().swap( m_hit_z     );
  std::vector<float>().swap( m_dca       );
  std::vector<float>().swap( m_hit_dE    );
  std::vector<float>().swap( m_hit_dE_dx );

}

StatusCode DCHDigiAlg::execute()
{
  info() << "Processing " << _nEvt << " events " << endmsg;
  Reset();
  std::map<unsigned long long, std::vector<edm4hep::SimTrackerHit> > id_hits_map;
  edm4hep::TrackerHitCollection* Vec   = w_DigiDCHCol.createAndPut();
  edm4hep::MCRecoTrackerAssociationCollection* AssoVec   = w_AssociationCol.createAndPut();
  const edm4hep::SimTrackerHitCollection* SimHitCol =  r_SimDCHCol.get();
  std::cout<<"input sim hit size="<< SimHitCol->size() <<std::endl;
  for( int i = 0; i < SimHitCol->size(); i++ ) 
  {
      edm4hep::SimTrackerHit SimHit = SimHitCol->at(i);
      unsigned long long id = SimHit.getCellID();
      float sim_hit_mom = sqrt( SimHit.getMomentum()[0]*SimHit.getMomentum()[0] + SimHit.getMomentum()[1]*SimHit.getMomentum()[1] + SimHit.getMomentum()[2]*SimHit.getMomentum()[2] );//GeV
      if(sim_hit_mom < m_mom_threshold) continue; 
      if(SimHit.getEDep() <= 0) continue; 
      
      if ( id_hits_map.find(id) != id_hits_map.end()) id_hits_map[id].push_back(SimHit);
      else 
      {
          std::vector<edm4hep::SimTrackerHit> vhit;
          vhit.push_back(SimHit);
          id_hits_map[id] = vhit ;
      }
  }

  for(std::map<unsigned long long, std::vector<edm4hep::SimTrackerHit> >::iterator iter = id_hits_map.begin(); iter != id_hits_map.end(); iter++)
  {
    unsigned long long wcellid = iter->first;
    auto trkHit = Vec->create();
    trkHit.setCellID(wcellid);
    double tot_edep   = 0 ;
    double tot_length = 0 ;
    double pos_x = 0 ;
    double pos_y = 0 ;
    double pos_z = 0 ;
    int simhit_size = iter->second.size();
    for(unsigned int i=0; i< simhit_size; i++)
    {
        tot_edep += iter->second.at(i).getEDep();//GeV
    }
    int chamber = m_decoder->get(wcellid, "chamber");
    int layer   = m_decoder->get(wcellid, "layer"  );
    int cellID  = m_decoder->get(wcellid, "cellID" );
    TVector3 Wstart(0,0,0);
    TVector3 Wend  (0,0,0);
    m_segmentation->cellposition(wcellid, Wstart, Wend);
    Wstart = 10*Wstart;// from DD4HEP cm to mm
    Wend   = 10*Wend  ;
    //std::cout<<"wcellid="<<wcellid<<",chamber="<<chamber<<",layer="<<layer<<",cellID="<<cellID<<",s_x="<<Wstart.x()<<",s_y="<<Wstart.y()<<",s_z="<<Wstart.z()<<",E_x="<<Wend.x()<<",E_y="<<Wend.y()<<",E_z="<<Wend.z()<<std::endl;

    TVector3  denominator = (Wend-Wstart) ;
    float min_distance = 999 ;
    for(unsigned int i=0; i< simhit_size; i++)
    {
        float sim_hit_mom = sqrt( iter->second.at(i).getMomentum()[0]*iter->second.at(i).getMomentum()[0] + iter->second.at(i).getMomentum()[1]*iter->second.at(i).getMomentum()[1] + iter->second.at(i).getMomentum()[2]*iter->second.at(i).getMomentum()[2] );//GeV
        float sim_hit_pt = sqrt( iter->second.at(i).getMomentum()[0]*iter->second.at(i).getMomentum()[0] + iter->second.at(i).getMomentum()[1]*iter->second.at(i).getMomentum()[1] );//GeV
        TVector3  pos(iter->second.at(i).getPosition()[0], iter->second.at(i).getPosition()[1], iter->second.at(i).getPosition()[2]);
        TVector3  numerator = denominator.Cross(Wstart-pos) ;
        float tmp_distance = numerator.Mag()/denominator.Mag() ;
        //std::cout<<"tmp_distance="<<tmp_distance<<",x="<<pos.x()<<",y="<<pos.y()<<",z="<<pos.z()<<",mom="<<sim_hit_mom<<",pt="<<sim_hit_pt<<std::endl;
        if(tmp_distance < min_distance){
            min_distance = tmp_distance;
            pos_x = pos.x();
            pos_y = pos.y();
            pos_z = pos.z();
        }
        tot_length += iter->second.at(i).getPathLength();//mm
 
        auto asso = AssoVec->create();
        asso.setRec(trkHit);
        asso.setSim(iter->second.at(i));
        asso.setWeight(iter->second.at(i).getEDep()/tot_edep);

        if(m_WriteAna){
            m_simhit_x.push_back(pos.x());
            m_simhit_y.push_back(pos.y());
            m_simhit_z.push_back(pos.z());
        }
    }
    
    trkHit.setTime(min_distance*1e3/m_velocity);//m_velocity is um/ns, drift time in ns
    trkHit.setEDep(tot_edep);// GeV
    trkHit.setEdx (tot_edep/tot_length); // GeV/mm
    trkHit.setPosition (edm4hep::Vector3d(pos_x, pos_y, pos_z));//position of closest sim hit
    trkHit.setCovMatrix(std::array<float, 6>{m_res_x, 0, m_res_y, 0, 0, m_res_z});//cov(x,x) , cov(y,x) , cov(y,y) , cov(z,x) , cov(z,y) , cov(z,z) in mm

    if(m_WriteAna){
        m_chamber.push_back(chamber);
        m_layer  .push_back(layer  );
        m_cell      .push_back(cellID);
        m_cell_x    .push_back(Wstart.x());
        m_cell_y    .push_back(Wstart.y());
        m_hit_x     .push_back(pos_x);
        m_hit_y     .push_back(pos_y);
        m_hit_z     .push_back(pos_z);
        m_dca       .push_back(min_distance);
        m_hit_dE    .push_back(trkHit.getEDep());
        m_hit_dE_dx .push_back(trkHit.getEdx() );
    }
  }
  std::cout<<"output digi DCHhit size="<< Vec->size() <<std::endl;
  _nEvt ++ ;

  if(m_WriteAna) m_tree->Fill();

  return StatusCode::SUCCESS;
}

StatusCode DCHDigiAlg::finalize()
{
  info() << "Processed " << _nEvt << " events " << endmsg;
  if(m_WriteAna){
      m_fout->cd();
      m_tree->Write();
      m_fout->Close();
  }
  return GaudiAlgorithm::finalize();
}
