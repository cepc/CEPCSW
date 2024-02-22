#ifndef Navigation_h
#define Navigation_h

#include "edm4hep/MCRecoTrackerAssociationCollection.h"
#include "edm4hep/TrackerHitCollection.h"
#include <map>

#if __has_include("edm4hep/EDM4hepVersion.h")
#include "edm4hep/EDM4hepVersion.h"
#else
// Copy the necessary parts from  the header above to make whatever we need to work here
#define EDM4HEP_VERSION(major, minor, patch) ((UINT64_C(major) << 32) | (UINT64_C(minor) << 16) | (UINT64_C(patch)))
// v00-09 is the last version without the capitalization change of the track vector members
#define EDM4HEP_BUILD_VERSION EDM4HEP_VERSION(0, 9, 0)
#endif

class Navigation{
 public:
  static Navigation* Instance();

  Navigation();
  ~Navigation();
  
  void Initialize();
  //void AddDataHandle(DataHandle* hdl){if(hdl)m_hdlVec.push_back(hdl);};
  void AddTrackerHitCollection(const edm4hep::TrackerHitCollection* col){m_hitColVec.push_back(col);};
  void AddTrackerAssociationCollection(const edm4hep::MCRecoTrackerAssociationCollection* col){m_assColVec.push_back(col);};

#if EDM4HEP_BUILD_VERSION <= EDM4HEP_VERSION(0, 10, 4)
  edm4hep::TrackerHit GetTrackerHit(const edm4hep::ObjectID& id, bool delete_by_caller=true);
  std::vector<edm4hep::SimTrackerHit> GetRelatedTrackerHit(const edm4hep::ObjectID& id);
#else
  edm4hep::TrackerHit GetTrackerHit(const podio::ObjectID& id, bool delete_by_caller=true);
  std::vector<edm4hep::SimTrackerHit> GetRelatedTrackerHit(const podio::ObjectID& id);
#endif
  std::vector<edm4hep::SimTrackerHit> GetRelatedTrackerHit(const edm4hep::TrackerHit& hit);
  std::vector<edm4hep::SimTrackerHit> GetRelatedTrackerHit(const edm4hep::TrackerHit& hit, const edm4hep::MCRecoTrackerAssociationCollection* col);
  
  //static Navigation* m_fNavigation;
 private:
  static Navigation* m_fNavigation;
  //DataHandle<edm4hep::MCRecoTrackerAssociationCollection> _inHitAssColHdl{"FTDStripTrackerHitsAssociation", Gaudi::DataHandle::Reader, this};
  std::vector<const edm4hep::TrackerHitCollection*> m_hitColVec;
  std::vector<const edm4hep::MCRecoTrackerAssociationCollection*> m_assColVec;
  std::map<int, edm4hep::TrackerHit> m_trkHits;
};
#endif 
