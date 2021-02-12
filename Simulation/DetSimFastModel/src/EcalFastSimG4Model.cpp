#include "EcalFastSimG4Model.h"

#include "G4Track.hh"
#include "G4FastTrack.hh"

EcalFastSimG4Model::EcalFastSimG4Model(G4String aModelName, G4Region* aEnvelope)
    : G4VFastSimulationModel(aModelName, aEnvelope) {

}

EcalFastSimG4Model::~EcalFastSimG4Model() {

}

G4bool EcalFastSimG4Model::IsApplicable(const G4ParticleDefinition& aParticle) {
//    return aParticle.GetPDGCharge() != 0;
    return true;
}

G4bool EcalFastSimG4Model::ModelTrigger(const G4FastTrack& aFastTrack) {
     //G4cout << __FILE__ << __LINE__ << ": ModelTrigger." << G4endl;

//    bool istrigged = false;
    bool istrigged = true;
    // only select the secondaries
    const G4Track* track = aFastTrack.GetPrimaryTrack();
    // secondaries
//G4cout << "trackID = " << track->GetTrackID() <<G4endl;
//    if (track->GetTrackID() != 0) {
//        istrigged = true;
//    }
//G4cout << "istrigged = " << istrigged <<G4endl;
    return istrigged;
}

void EcalFastSimG4Model::DoIt(const G4FastTrack& aFastTrack, G4FastStep& aFastStep) {
     //G4cout << __FILE__ << __LINE__ << ": DoIt." << G4endl;

    aFastStep.ProposeTrackStatus(fStopAndKill);
}
