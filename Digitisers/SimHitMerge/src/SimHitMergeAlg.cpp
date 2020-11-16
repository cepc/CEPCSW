#include "SimHitMergeAlg.h"

#include "edm4hep/SimCalorimeterHit.h"
#include "edm4hep/Vector3f.h"

#include "DD4hep/Detector.h"
#include "DD4hep/IDDescriptor.h"
#include "DD4hep/Plugins.h"

#include <string>
#include <iostream>
#include <cmath>
#include <stdexcept>
#include <sstream>



DECLARE_COMPONENT( SimHitMergeAlg )






SimHitMergeAlg::SimHitMergeAlg(const std::string& name, ISvcLocator* svcLoc)
     : GaudiAlgorithm(name, svcLoc)
      
{
}

StatusCode SimHitMergeAlg::initialize() {
     m_geosvc = service<IGeomSvc>("GeomSvc");
     if (!m_geosvc) {
       error() << "Failed to find GeomSvc." << endmsg;
       return StatusCode::FAILURE;
     }
     m_dd4hep_geo = m_geosvc->lcdd();
     if (!m_dd4hep_geo) {
       error() << "failed to retrieve dd4hep_geo: " << m_dd4hep_geo << endmsg;
       return StatusCode::FAILURE;
     }
     m_cellIDConverter = new dd4hep::rec::CellIDPositionConverter(*m_dd4hep_geo);

     if(m_inputColNames.size() != m_outputColNames.size() ) throw ("Error input size != output size");
     for (auto& input : m_inputColNames) {
	  m_InputCollections.push_back( new SimCaloType(input, Gaudi::DataHandle::Reader, this) );
     }

     for (auto& output : m_outputColNames) {
	  m_OutputCollections.push_back( new SimCaloType(output, Gaudi::DataHandle::Writer, this) );
     }
     

     return GaudiAlgorithm::initialize();
}

StatusCode SimHitMergeAlg::execute()
{

     for (unsigned int k0 = 0; k0 < m_inputColNames.size(); k0++)
     {
          std::map<unsigned long long, edm4hep::SimCalorimeterHit> id_hit_map;
          std::map<unsigned long long, std::vector<edm4hep::ConstCaloHitContribution> > id_vconb_map;
	  edm4hep::SimCalorimeterHitCollection* mergedCol = m_OutputCollections[k0]->createAndPut();
	  auto col = m_InputCollections[k0]->get();
          //std::cout<<"input="<<m_InputCollections[k0]->objKey()<<",size="<<col->size()<<std::endl;
	  for (auto Simhit: *col){
              auto id = Simhit.getCellID();
              if(Simhit.getEnergy() <=0 ) continue; 
              if ( id_hit_map.find(id) != id_hit_map.end()) id_hit_map[id].setEnergy(id_hit_map[id].getEnergy() + Simhit.getEnergy());
              else id_hit_map[id] = Simhit;
              std::vector<edm4hep::ConstCaloHitContribution> tmp_vconb ;
              for(int kk=0; kk<Simhit.contributions_size(); kk++){
                  tmp_vconb.push_back(Simhit.getContributions(kk));
              }
              
              if( id_vconb_map.find(id) != id_vconb_map.end()){
                  for(int kk=0; kk<tmp_vconb.size();kk++){
                      id_vconb_map[id].push_back(tmp_vconb.at(kk));
                  }
              }
              else id_vconb_map[id] = tmp_vconb;
                  
          }
        
          for(std::map<unsigned long long, edm4hep::SimCalorimeterHit>::iterator iter = id_hit_map.begin(); iter != id_hit_map.end(); iter++)
          {
	       edm4hep::SimCalorimeterHit Simhit = iter->second ;
               dd4hep::Position position = m_cellIDConverter->position(Simhit.getCellID());//cm
	       edm4hep::Vector3f hitPos(position.x()*10, position.y()*10, position.z()*10);//to mm
               //std::cout<<"id="<<Simhit.getCellID()<<",hitPos.x="<<hitPos[0]<<",y="<<hitPos[1]<<",z="<<hitPos[2]<<",ori x="<<Simhit.getPosition()[0]<<",y="<<Simhit.getPosition()[1]<<",z="<<Simhit.getPosition()[2]<<std::endl;
	       auto Mergedhit = mergedCol->create();
	       Mergedhit.setCellID  (Simhit.getCellID());
	       Mergedhit.setPosition(hitPos);
	       Mergedhit.setEnergy  (Simhit.getEnergy());
               for(int ii=0; ii<id_vconb_map[iter->first].size(); ii++){
	           Mergedhit.addToContributions(id_vconb_map[iter->first].at(ii));
               }
          }
          /*
          std::cout<<"output="<<m_OutputCollections[k0]->objKey()<<",size="<<mergedCol->size()<<std::endl;
          for(int ii=0; ii<mergedCol->size(); ii++){
              auto tmp_hit = mergedCol->at(ii); 
              std::cout<<"id="<<tmp_hit.getCellID()<<",E="<<tmp_hit.getEnergy()<<",conb size="<<tmp_hit.contributions_size()<<std::endl;
          }
          */
          
     }     
     return StatusCode::SUCCESS;
}

StatusCode SimHitMergeAlg::finalize()
{
     std::cout<<"SimHitMergeAlg FINISHED"<<std::endl;


     return GaudiAlgorithm::finalize();
}
