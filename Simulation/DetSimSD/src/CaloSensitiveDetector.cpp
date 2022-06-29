#include "DetSimSD/CaloSensitiveDetector.h"

#include "G4SDManager.hh"

#include <algorithm>

CaloSensitiveDetector::CaloSensitiveDetector(const std::string& name, dd4hep::Detector& description, bool is_merge_enabled)
    : DDG4SensitiveDetector(name, description),
      m_hc(nullptr),
      m_isMergeEnabled(is_merge_enabled){
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
    if(m_applyBirksLaw) h.doApplyBirksLaw();
    dd4hep::Position pos = 0.5 * (h.prePos() + h.postPos());
    HitContribution contrib = dd4hep::sim::Geant4Hit::extractContribution(step);
    const std::string& name = GetName();
    unsigned long id = getCellID( step );
    CalorimeterHit* hit=nullptr;
    if(m_isMergeEnabled){
      std::map<unsigned long, CalorimeterHit*>::iterator it = m_hitMap.find(id);
      if(it!=m_hitMap.end()) hit = it->second;
    }
    else{
      //Commented by fucd: hit position almost different, only very few hits found sucessfully, so discard to find since this option is disable merge
      //hit=find(m_hc,dd4hep::sim::HitPositionCompare<CalorimeterHit>(pos));
    }
    //    G4cout << "----------- Geant4GenericSD<Calorimeter>::buildHits : position : " << pos << G4endl;
    if ( !hit ) {
        // not applicable for segmentation case
        //G4ThreeVector local(0,0,0);
        //G4ThreeVector global = h.preTouchable()->GetHistory()->GetTopTransform().InverseTransformPoint(local);
        //hit = new CalorimeterHit(dd4hep::Position(global.x(), global.y(), global.z()));
        if(m_isMergeEnabled){
          dd4hep::Position posCellCenter = getNominalPosition(step, id);
          hit = new CalorimeterHit(posCellCenter);
          m_hitMap[id] = hit;
        }
        else hit = new CalorimeterHit(pos);
        hit->cellID  = id; //getCellID( step );
        m_hc->insert(hit);
    }
    hit->truth.push_back(contrib);
    //hit->energyDeposit += contrib.deposit;
    hit->energyDeposit += h.totalEnergy();
    //std::cout << "Apply Birk law: before = " << contrib.deposit << " after = " << h.totalEnergy() << std::endl;
    
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
