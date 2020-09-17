#include "DataHelper/Navigation.h"

#include "edm4hep/SimTrackerHit.h"
#include "edm4hep/TrackerHit.h"

Navigation* Navigation::m_fNavigation = nullptr;

Navigation* Navigation::Instance(){
  if(!m_fNavigation) m_fNavigation = new Navigation();
  return m_fNavigation;
}

Navigation::Navigation(){
}

Navigation::~Navigation(){
}

void Navigation::Initialize(){
  m_hitColVec.clear();
  m_assColVec.clear();
  for(std::map<int, edm4hep::ConstTrackerHit>::iterator it=m_trkHits.begin();it!=m_trkHits.end();it++){
    // delete it->second;
  }
  m_trkHits.clear();
}

edm4hep::ConstTrackerHit Navigation::GetTrackerHit(const edm4hep::ObjectID& obj_id, bool delete_by_caller){
  int id = obj_id.collectionID * 10000000 + obj_id.index;
  if(!delete_by_caller){
    if(m_trkHits.find(id)!=m_trkHits.end()) return m_trkHits[id];
  }
  /*
  for(int i=0;i<m_assColVec.size();i++){
    for(auto ass : *m_assColVec[i]){
      edm4hep::ObjectID rec_id = ass.getRec().getObjectID();
      if(rec_id.collectionID!=id.collectionID)break;
      else if(rec_id.index==id.index){
	m_trkHits.push_back(ass.getRec());
	return &(m_trkHits.back());
      }
    }
  }
  */
  for(int i=0;i<m_hitColVec.size();i++){
    for(auto hit : *m_hitColVec[i]){
      edm4hep::ObjectID this_id = hit.getObjectID();
      if(this_id.collectionID!=obj_id.collectionID)break;
      else if(this_id.index==obj_id.index){
	edm4hep::ConstTrackerHit hit_copy = edm4hep::ConstTrackerHit(hit);
	if(!delete_by_caller) m_trkHits[id] = hit_copy;
	return hit_copy;//&(m_trkHits[id]);
      }
    }
  }
  
  throw std::runtime_error("Not found TrackerHit");
}

std::vector<edm4hep::ConstSimTrackerHit> Navigation::GetRelatedTrackerHit(const edm4hep::ObjectID& id){
  std::vector<edm4hep::ConstSimTrackerHit> hits;
  for(int i=0;i<m_assColVec.size();i++){
    for(auto ass : *m_assColVec[i]){
      edm4hep::ObjectID this_id = ass.getRec().getObjectID();
      if(this_id.collectionID!=id.collectionID)break;
      else if(this_id.index==id.index) hits.push_back(ass.getSim()); 
    }
  }
  return hits;
}

std::vector<edm4hep::ConstSimTrackerHit> Navigation::GetRelatedTrackerHit(const edm4hep::TrackerHit& hit){
  std::vector<edm4hep::ConstSimTrackerHit> hits;
  for(int i=0;i<m_assColVec.size();i++){
    for(auto ass : *m_assColVec[i]){
      if(ass.getRec().getObjectID().collectionID != hit.getObjectID().collectionID) break;
      else if(ass.getRec()==hit) hits.push_back(ass.getSim());
    }
  }
  return hits;
}
