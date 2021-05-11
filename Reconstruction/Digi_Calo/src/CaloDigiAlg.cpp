/* -*- Mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */
#include "CaloDigiAlg.h"


#include "edm4hep/SimCalorimeterHit.h"
#include "edm4hep/CalorimeterHit.h"
#include "edm4hep/Vector3f.h"

#include "DD4hep/Detector.h"
#include <DD4hep/Objects.h>


#include <math.h>
#include <cmath>
#include <algorithm>

DECLARE_COMPONENT( CaloDigiAlg )

CaloDigiAlg::CaloDigiAlg(const std::string& name, ISvcLocator* svcLoc)
  : GaudiAlgorithm(name, svcLoc),
    _nEvt(0)
{
  
  // Input collections
  declareProperty("SimCaloHitCollection", r_SimCaloCol, "Handle of the Input SimCaloHit collection");
  
  // Output collections
  declareProperty("CaloHitCollection", w_DigiCaloCol, "Handle of Digi CaloHit collection");
  
  declareProperty("CaloAssociationCollection", w_CaloAssociationCol, "Handle of CaloAssociation collection");
   
}

StatusCode CaloDigiAlg::initialize()
{

  std::cout<<"CaloDigiAlg::m_scale="<<m_scale<<std::endl;
  m_geosvc = service<IGeomSvc>("GeomSvc");
  if ( !m_geosvc )  throw "CaloDigiAlg :Failed to find GeomSvc ...";
  dd4hep::Detector* m_dd4hep = m_geosvc->lcdd();
  if ( !m_dd4hep )  throw "CaloDigiAlg :Failed to get dd4hep::Detector ...";
  m_cellIDConverter = new dd4hep::rec::CellIDPositionConverter(*m_dd4hep);
  return GaudiAlgorithm::initialize();
}

StatusCode CaloDigiAlg::execute()
{
  std::map<unsigned long long, edm4hep::SimCalorimeterHit> id_hit_map;
  std::map<unsigned long long, std::vector<edm4hep::SimCalorimeterHit> > id_hits_map;
  edm4hep::CalorimeterHitCollection* caloVec   = w_DigiCaloCol.createAndPut();
  edm4hep::MCRecoCaloAssociationCollection* caloAssoVec   = w_CaloAssociationCol.createAndPut();
  const edm4hep::SimCalorimeterHitCollection* SimHitCol =  r_SimCaloCol.get();
  double tot_e = 0 ;
  if(SimHitCol == 0) 
  {
     std::cout<<"not found SimCalorimeterHitCollection"<< std::endl;
     return StatusCode::SUCCESS;
  }
  std::cout<<"digi, input sim hit size="<< SimHitCol->size() <<std::endl;
  for( int i = 0; i < SimHitCol->size(); i++ ) 
  {
      edm4hep::SimCalorimeterHit SimHit = SimHitCol->at(i);
      unsigned long long id = SimHit.getCellID();
      float en = SimHit.getEnergy();
      tot_e += en;
      if ( id_hit_map.find(id) != id_hit_map.end()) id_hit_map[id].setEnergy(id_hit_map[id].getEnergy() + en);
      else id_hit_map[id] = SimHit ;

      if ( id_hits_map.find(id) != id_hits_map.end()) id_hits_map[id].push_back(SimHit);
      else 
      {
          std::vector<edm4hep::SimCalorimeterHit> vhit;
          vhit.push_back(SimHit);
          id_hits_map[id] = vhit ;
      }
  }
  for(std::map<unsigned long long, edm4hep::SimCalorimeterHit>::iterator iter = id_hit_map.begin(); iter != id_hit_map.end(); iter++)
  {
    auto caloHit = caloVec->create();
    caloHit.setCellID((iter->second).getCellID());
    caloHit.setEnergy((iter->second).getEnergy()*m_scale);
    dd4hep::Position position = m_cellIDConverter->position(caloHit.getCellID());
    edm4hep::Vector3f vpos(position.x()*10, position.y()*10, position.z()*10);// cm to mm
    caloHit.setPosition(vpos);
    //std::cout << "sim hit id =" << caloHit.getCellID() <<",x="<<position.x()<<",y="<<position.y()<<",z="<<position.z() <<",real x="<<(iter->second).getPosition().x <<",y="<<(iter->second).getPosition().y<<",z="<<(iter->second).getPosition().z<< std::endl;
    
    if( id_hits_map.find(iter->first) != id_hits_map.end()) 
    {
        for(unsigned int i=0; i< id_hits_map[iter->first].size(); i++)
        {
            auto asso = caloAssoVec->create();
            asso.setRec(caloHit);
            asso.setSim(id_hits_map[iter->first].at(i));
            asso.setWeight(id_hits_map[iter->first].at(i).getEnergy()/(iter->second).getEnergy());
        }
    }
    else std::cout<<"Error in Digi Calo"<<std::endl;
  }
    
  std::cout<<"total sim e ="<< tot_e <<std::endl;
  std::cout<<"digi, output digi hit size="<< caloVec->size() <<std::endl;
  std::cout<<"digi, output caloAssoVec hit size="<< caloAssoVec->size() <<std::endl;
  _nEvt ++ ;

  return StatusCode::SUCCESS;
}

StatusCode CaloDigiAlg::finalize()
{
  info() << "Processed " << _nEvt << " events " << endmsg;
  return GaudiAlgorithm::finalize();
}
