/* -*- Mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */
#include "DCHDigiAlg.h"


#include "edm4hep/SimCalorimeterHit.h"
#include "edm4hep/CalorimeterHit.h"
#include "edm4hep/Vector3f.h"

#include "DD4hep/Detector.h"
#include <DD4hep/Objects.h>


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
  
  declareProperty("CaloAssociationCollection", w_AssociationCol, "Handle of Association collection");
   
}

StatusCode DCHDigiAlg::initialize()
{
  /*
  m_geosvc = service<IGeoSvc>("GeoSvc");
  if ( !m_geosvc )  throw "DCHDigiAlg :Failed to find GeoSvc ...";
  dd4hep::Detector* m_dd4hep = m_geosvc->lcdd();
  if ( !m_dd4hep )  throw "DCHDigiAlg :Failed to get dd4hep::Detector ...";
  m_cellIDConverter = new dd4hep::rec::CellIDPositionConverter(*m_dd4hep);
  */
  std::cout<<"DCHDigiAlg::initialized"<< std::endl;
  return GaudiAlgorithm::initialize();
}

StatusCode DCHDigiAlg::execute()
{
  std::map<unsigned long long, std::vector<edm4hep::SimTrackerHit> > id_hits_map;
  edm4hep::TrackerHitCollection* Vec   = w_DigiDCHCol.createAndPut();
  edm4hep::MCRecoTrackerAssociationCollection* AssoVec   = w_AssociationCol.createAndPut();
  const edm4hep::SimTrackerHitCollection* SimHitCol =  r_SimDCHCol.get();
  if(SimHitCol == 0) 
  {
     std::cout<<"not found SimCalorimeterHitCollection"<< std::endl;
     return StatusCode::SUCCESS;
  }
  std::cout<<"input sim hit size="<< SimHitCol->size() <<std::endl;
  for( int i = 0; i < SimHitCol->size(); i++ ) 
  {
      edm4hep::SimTrackerHit SimHit = SimHitCol->at(i);
      unsigned long long id = SimHit.getCellID();
      
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
    auto trkHit = Vec->create();
    trkHit.setCellID(iter->first);
    double tot_edep   = 0 ;
    double tot_length = 0 ;
    double tot_time = 0 ;
    double tot_x = 0 ;
    double tot_y = 0 ;
    double tot_z = 0 ;
    int simhit_size = iter->second.size();
    for(unsigned int i=0; i< simhit_size; i++)
    {
        tot_edep   += iter->second.at(i).getEDep();//GeV
        tot_length += iter->second.at(i).getPathLength();//mm
        tot_time += iter->second.at(i).getTime();
        tot_x    += iter->second.at(i).getPosition()[0];
        tot_y    += iter->second.at(i).getPosition()[1];
        tot_z    += iter->second.at(i).getPosition()[2];
 
        auto asso = AssoVec->create();
        asso.setRec(trkHit);
        asso.setSim(iter->second.at(i));
        asso.setWeight(1.0/simhit_size);
    }
    
    trkHit.setTime(tot_time/simhit_size);
    trkHit.setEDep(tot_edep/simhit_size);
    trkHit.setEdx (tot_edep*1000/(tot_length/10) ); // MeV/cm
    trkHit.setPosition (edm4hep::Vector3d(tot_x/simhit_size, tot_y/simhit_size, tot_z/simhit_size));
  }
  std::cout<<"output digi DCHhit size="<< Vec->size() <<std::endl;
  _nEvt ++ ;

  return StatusCode::SUCCESS;
}

StatusCode DCHDigiAlg::finalize()
{
  info() << "Processed " << _nEvt << " events " << endmsg;
  return GaudiAlgorithm::finalize();
}
