#include "DetSimSD/CaloSensitiveDetector.h"

#include "G4SDManager.hh"

#include <algorithm>

CaloSensitiveDetector::CaloSensitiveDetector(const std::string& name, dd4hep::Detector& description)
    : DDG4SensitiveDetector(name, description),
      m_hc(nullptr) {
    const std::string& coll_name = m_sensitive.hitsCollection();
    collectionName.insert(coll_name);
}

void
CaloSensitiveDetector::Initialize(G4HCofThisEvent* HCE) {

    // the collection name is provided by DD4hep
    const std::string& coll_name = collectionName[0];
    // m_hc = new G4THitsCollection<CalorimeterHit>(GetName(), coll_name);
    m_hc = new HitCollection(GetName(), coll_name);

    int HCID = -1;
    if(HCID<0) HCID = G4SDManager::GetSDMpointer()->GetCollectionID(m_hc);
    HCE->AddHitsCollection( HCID, m_hc ); 

    m_hitMap.clear();
}

G4bool
CaloSensitiveDetector::ProcessHits(G4Step* step, G4TouchableHistory*) {

    // std::cout << "CaloSensitiveDetector::ProcessHits" << std::endl;

    dd4hep::sim::Geant4StepHandler h(step);
    dd4hep::Position pos = 0.5 * (h.prePos() + h.postPos());
    HitContribution contrib = dd4hep::sim::Geant4Hit::extractContribution(step);
    const std::string& name = GetName();
    unsigned long id = getCellID( step );
    //std::cout << name << " " << id << std::endl;
    CalorimeterHit* hit=nullptr;
    if(name.find("Merge")==-1&&name.find("merge")==-1) find(m_hc,dd4hep::sim::HitPositionCompare<CalorimeterHit>(pos));
    else{
      std::map<unsigned long, CalorimeterHit*>::iterator it = m_hitMap.find(id);
      if(it!=m_hitMap.end()) hit = it->second;
    }
    //    G4cout << "----------- Geant4GenericSD<Calorimeter>::buildHits : position : " << pos << G4endl;
    if ( !hit ) {
        hit = new CalorimeterHit(pos);
        hit->cellID  = getCellID( step );
        m_hc->insert(hit);
	if(name.find("Merge")!=-1||name.find("merge")!=-1) m_hitMap[id] = hit;
    }
    hit->truth.push_back(contrib);
    hit->energyDeposit += contrib.deposit;


    
    return true;
}

void
CaloSensitiveDetector::EndOfEvent(G4HCofThisEvent* HCE) {

}

CaloSensitiveDetector::CalorimeterHit*
CaloSensitiveDetector::find(const CaloSensitiveDetector::HitCollection* c,
                            const dd4hep::sim::HitCompare<CaloSensitiveDetector::CalorimeterHit>& cmp) {
    typedef std::vector<CalorimeterHit*> _V;
    const _V* v = (const _V*) c->GetVector();
    for (_V::const_iterator i = v->begin(); i != v->end(); ++i) {
        if (cmp(*i)) {
            return *i;
        }
    }
    return 0;
}
