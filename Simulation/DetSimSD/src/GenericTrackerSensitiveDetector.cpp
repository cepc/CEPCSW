#include "GenericTrackerSensitiveDetector.h"

#include "G4Step.hh"
#include "G4VProcess.hh"
//#include "G4SteppingManager.hh"
#include "G4SDManager.hh"
//#include "UserTrackInformation.hh"
#include "DD4hep/DD4hepUnits.h"

GenericTrackerSensitiveDetector::GenericTrackerSensitiveDetector(const std::string& name,
								 dd4hep::Detector& description)
  : DDG4SensitiveDetector(name, description),
    m_hc(nullptr){
  
  G4String CollName=name+"Collection";
  collectionName.insert(CollName);
}

void GenericTrackerSensitiveDetector::Initialize(G4HCofThisEvent* HCE){
  m_hc = new HitCollection(GetName(), collectionName[0]);
  int HCID = G4SDManager::GetSDMpointer()->GetCollectionID(m_hc);
  HCE->AddHitsCollection( HCID, m_hc ); 
}

G4bool GenericTrackerSensitiveDetector::ProcessHits(G4Step* step, G4TouchableHistory*){
  
  G4TouchableHandle touchPost = step->GetPostStepPoint()->GetTouchableHandle(); 
  G4TouchableHandle touchPre  = step->GetPreStepPoint()->GetTouchableHandle(); 
  dd4hep::sim::Geant4StepHandler h(step);
  dd4hep::Position prePos    = h.prePos();
  dd4hep::Position postPos   = h.postPos();
  dd4hep::Position direction = postPos - prePos;
  dd4hep::Position position  = mean_direction(prePos,postPos);
  double   hit_len   = direction.R();
  if (hit_len > 0) {
    double new_len = mean_length(h.preMom(),h.postMom())/hit_len;
    direction *= new_len/hit_len;
  }
  dd4hep::sim::Geant4TrackerHit* hit = nullptr;
  hit = new dd4hep::sim::Geant4TrackerHit(h.track->GetTrackID(),
					  h.track->GetDefinition()->GetPDGEncoding(),
					  step->GetTotalEnergyDeposit(),
					  h.track->GetGlobalTime());

  if ( hit )  {
    hit->cellID  = getCellID( step ) ;
    hit->energyDeposit =  step->GetTotalEnergyDeposit();
    hit->position = position;
    hit->momentum = direction;
    hit->length   = hit_len;
    m_hc->insert(hit);
    return true;
  }
  throw std::runtime_error("new() failed: Cannot allocate hit object");
  return false;
}

void GenericTrackerSensitiveDetector::EndOfEvent(G4HCofThisEvent* HCE){
}
