#include "TrackerCombineSensitiveDetector.h"

#include "G4Step.hh"
#include "G4VProcess.hh"
#include "G4SDManager.hh"
#include "DD4hep/DD4hepUnits.h"

TrackerCombineSensitiveDetector::TrackerCombineSensitiveDetector(const std::string& name,
								 dd4hep::Detector& description)
  : DDG4SensitiveDetector(name, description),
    m_hc(nullptr){
  G4String CollName = m_sensitive.hitsCollection();
  collectionName.insert(CollName);
}

void TrackerCombineSensitiveDetector::Initialize(G4HCofThisEvent* HCE){
  userData.e_cut = m_sensitive.energyCutoff();

  m_hc = new HitCollection(GetName(), collectionName[0]);
  int HCID = G4SDManager::GetSDMpointer()->GetCollectionID(m_hc);
  HCE->AddHitsCollection( HCID, m_hc ); 
}

G4bool TrackerCombineSensitiveDetector::ProcessHits(G4Step* step, G4TouchableHistory*){
  dd4hep::sim::Geant4StepHandler h(step);
  bool return_code = false;
  if ( userData.current == -1 ) userData.start(getCellID(step), step, h.pre);
  else if ( !userData.track || userData.current != h.track->GetTrackID() ) {
    return_code = userData.extractHit(m_hc) != 0;
    userData.start(getCellID(step), step, h.pre);
  }

  // ....update .....
  userData.update(step);

  void *prePV = h.volume(h.pre), *postPV = h.volume(h.post);
  if ( prePV != postPV ) {
    return_code = userData.extractHit(m_hc) != 0;
    void* postSD = h.sd(h.post);
    if ( 0 != postSD )   {
      void* preSD = h.sd(h.pre);
      if ( preSD == postSD ) {
	// fucd: getCellID(step) for preVolume not postVolume, so should start at next step 
	//userData.start(getCellID(step), step, h.post);
      }
    }
  }
  else if ( userData.track->GetTrackStatus() == fStopAndKill ) {
    return_code = userData.extractHit(m_hc) != 0;
  }
  return return_code;
}

void TrackerCombineSensitiveDetector::EndOfEvent(G4HCofThisEvent* HCE){
  userData.clear();
}
