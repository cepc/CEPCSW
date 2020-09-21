#include "DriftChamberSensitiveDetector.h"

#include "G4SDManager.hh"

DriftChamberSensitiveDetector::DriftChamberSensitiveDetector(const std::string& name,
                                                             dd4hep::Detector& description)
    : DDG4SensitiveDetector(name, description),
      m_hc(nullptr) {

    const std::string& coll_name = m_sensitive.hitsCollection();

    collectionName.insert(coll_name);
}

bool DriftChamberSensitiveDetector::setDedxSimTool(ToolHandle<IDedxSimTool> simtool) {
    m_dedx_simtool = simtool;

    return true;
}

void
DriftChamberSensitiveDetector::Initialize(G4HCofThisEvent* HCE) {

    const std::string& coll_name = collectionName[0];
    m_hc = new HitCollection(GetName(), coll_name);

    int HCID = -1;
    if(HCID<0) HCID = G4SDManager::GetSDMpointer()->GetCollectionID(m_hc);
    HCE->AddHitsCollection( HCID, m_hc ); 

}

G4bool
DriftChamberSensitiveDetector::ProcessHits(G4Step* step, G4TouchableHistory*) {
    // Refer to: DDG4/legacy/Geant4TrackerSD.cpp (note: there's bug in momentum calculation)
    //           DDCore/include/DD4hep/Objects.h (mean_direction and mean_length)

    dd4hep::sim::Geant4StepHandler h(step);

    dd4hep::Position prePos    = h.prePos();
    dd4hep::Position postPos   = h.postPos();
    dd4hep::Position direction = postPos - prePos;
    dd4hep::Position position  = mean_direction(prePos,postPos); // (pre+post)/2
    double           hit_len   = direction.R();

    HitContribution contrib = dd4hep::sim::Geant4Hit::extractContribution(step);
    // Now, invokes the dE/dx simulator
    double dedx = 0.0;
    dedx = m_dedx_simtool->dedx(step);
    // G4cout << "-----> dedx: " << dedx << G4endl;

    double de = hit_len * dedx;
    // contrib.deposit = de; // if need the de from dedx simulator

    // create a new hit
    TrackerHit* hit = new TrackerHit(
                                     h.track->GetTrackID(),
                                     h.track->GetDefinition()->GetPDGEncoding(),
                                     de, // not the Geant4's deposit energy. from dE/dx simulator
                                     h.track->GetGlobalTime()
                                     );
    hit->cellID = getCellID(step);
    hit->energyDeposit = de; // FIXME: also use the dedx
    hit->position = position;
    hit->momentum = (h.preMom() + h.postMom() )/2;
    hit->length   = hit_len;
    m_hc->insert(hit);

    return true;
}

void
DriftChamberSensitiveDetector::EndOfEvent(G4HCofThisEvent* HCE) {

}
