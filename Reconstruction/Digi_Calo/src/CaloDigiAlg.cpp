/* -*- Mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */
#include "CaloDigiAlg.h"


#include "edm4hep/SimCalorimeterHit.h"
#include "edm4hep/CalorimeterHit.h"
#include "edm4hep/Vector3f.h"

#include "DD4hep/Detector.h"
#include <DD4hep/Objects.h>
#include <DDRec/DetectorData.h>
#include <DDRec/CellIDPositionConverter.h>


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
  
   
}

StatusCode CaloDigiAlg::initialize()
{

  std::cout<<"CaloDigiAlg::m_scale="<<m_scale<<std::endl;
  m_geosvc = service<IGeoSvc>("GeoSvc");
  if ( !m_geosvc )  throw "CaloDigiAlg :Failed to find GeoSvc ...";
  dd4hep::Detector* m_dd4hep = m_geosvc->lcdd();
  if ( !m_dd4hep )  throw "CaloDigiAlg :Failed to get dd4hep::Detector ...";
  m_cellIDConverter = new dd4hep::rec::CellIDPositionConverter(*m_dd4hep);
  try{
      const dd4hep::DetElement &detElement = m_dd4hep->detector("CaloDetector");
      dd4hep::rec::LayeredCalorimeterData* Data = detElement.extension<dd4hep::rec::LayeredCalorimeterData>() ;
      const std::vector<dd4hep::rec::LayeredCalorimeterData::Layer>& layerLayout = Data->layers;
      assert(layerLayout.size()>=1);
      m_layerLayout = &layerLayout;
      std::cout<<"saved m_layerLayout"<<std::endl;
  }
  catch(...) 
  {
      throw "CaloDigiAlg :Failed to get LayeredCalorimeterData ...";
  }
  return GaudiAlgorithm::initialize();
}

StatusCode CaloDigiAlg::execute()
{
  std::map<std::string, std::vector<edm4hep::SimCalorimeterHit> > id_vhit_map;
  edm4hep::CalorimeterHitCollection* caloVec   = w_DigiCaloCol.createAndPut();
  const edm4hep::SimCalorimeterHitCollection* SimHitCol =  r_SimCaloCol.get();
  if(SimHitCol == 0) 
  {
     std::cout<<"not found SimCalorimeterHitCollection"<< std::endl;
     return StatusCode::SUCCESS;
  }
  std::cout<<"digi, input sim hit size="<< SimHitCol->size() <<std::endl;
  float tot_sim_en = 0;
  for( int i = 0; i < SimHitCol->size(); i++ ) 
  {
      edm4hep::SimCalorimeterHit SimHit = SimHitCol->at(i);
      float en = SimHit.getEnergy();
      tot_sim_en += en;
      std::string id = GetCode(&SimHit);
      //std::cout<<"sim i ="<< i<<",id="<< id<<",x="<<SimHit.getPosition().x<<",y="<<SimHit.getPosition().y<<",z="<<SimHit.getPosition().z <<std::endl;
      if(id_vhit_map.find(id) != id_vhit_map.end()) id_vhit_map[id].push_back(SimHit);
      else
      {
          std::vector<edm4hep::SimCalorimeterHit> vect;
          vect.push_back(SimHit);
          id_vhit_map[id] = vect;
      }
  }
  float tot_en = 0;
  for(std::map<std::string, std::vector <edm4hep::SimCalorimeterHit> >::iterator iter = id_vhit_map.begin(); iter != id_vhit_map.end(); iter++)
  {
    float energy = 0;
    float x0 = 0;
    float y0 = 0;
    float z0 = 0;
    for(int j=0; j < iter->second.size(); j++)
    {
        float hit_en = iter->second.at(j).getEnergy();
        energy += hit_en;
        x0     += hit_en*iter->second.at(j).getPosition().x ;
        y0     += hit_en*iter->second.at(j).getPosition().y ;
        z0     += hit_en*iter->second.at(j).getPosition().z ;
    }
    x0 = x0/energy;
    y0 = y0/energy;
    z0 = z0/energy;
    edm4hep::Vector3f pos(x0, y0, z0);
    energy = energy*m_scale;
    auto caloHit = caloVec->create();
    //caloHit.setCellID  (iter->second.getCellID());
    caloHit.setEnergy  (energy);
    caloHit.setPosition(pos);
    tot_en += caloHit.getEnergy();
  }
    
  std::cout<<"total sim e ="<< tot_sim_en <<std::endl;
  std::cout<<"digi, output digi hit size="<< caloVec->size()<<",tot_en="<<tot_en <<std::endl;
  _nEvt ++ ;

  return StatusCode::SUCCESS;
}

StatusCode CaloDigiAlg::finalize()
{
  info() << "Processed " << _nEvt << " events " << endmsg;
  return GaudiAlgorithm::finalize();
}


std::string CaloDigiAlg::GetCode(const edm4hep::SimCalorimeterHit *const pHit) const
{
    int layer = -1 ;   
    for (unsigned int i = 0, iMax = m_layerLayout->size(); i < iMax; ++i)
    {
        const float distance(m_layerLayout->at(i).distance);
        const float sensitive_thickness(m_layerLayout->at(i).sensitive_thickness);
        if( (distance - 0.5*sensitive_thickness) <= pHit->getPosition().x && pHit->getPosition().x <= (distance + 0.5*sensitive_thickness)) {layer = i ; break;}
    }
    if(layer==-1)
    {
        int lmax = m_layerLayout->size()-1;
        std::cout<<"Error BarrelLayer, set to default 1, Hit.x="<<pHit->getPosition().x<<", min_x="<<m_layerLayout->at(0).distance-0.5*m_layerLayout->at(0).sensitive_thickness<<", max_x="<<m_layerLayout->at(lmax).distance+0.5*m_layerLayout->at(lmax).sensitive_thickness<<std::endl;
        layer = 0;
    }
    float cellSize0     = m_layerLayout->at(layer).cellSize0;
    float cellSize1     = m_layerLayout->at(layer).cellSize1;
    int cell0 = floor((pHit->getPosition().y+0.5*cellSize0)/cellSize0);
    int cell1 = floor((pHit->getPosition().z+0.5*cellSize1)/cellSize1);
    std::string sl = std::to_string(layer);
    std::string s0 = std::to_string(cell0);
    std::string s1 = std::to_string(cell1);
    std::string id = sl+"_"+s0+"_"+s1 ;
    return id;
}
