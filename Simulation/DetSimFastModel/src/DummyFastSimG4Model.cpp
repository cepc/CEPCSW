#include "DummyFastSimG4Model.h"

#include "G4Track.hh"
#include "G4FastTrack.hh"

DummyFastSimG4Model::DummyFastSimG4Model(G4String aModelName, G4Region* aEnvelope)
    : G4VFastSimulationModel(aModelName, aEnvelope) {

}

DummyFastSimG4Model::~DummyFastSimG4Model() {

}

G4bool DummyFastSimG4Model::IsApplicable(const G4ParticleDefinition& aParticle) {
    return aParticle.GetPDGCharge() != 0;
}

G4bool DummyFastSimG4Model::ModelTrigger(const G4FastTrack& aFastTrack) {
    // G4cout << __FILE__ << __LINE__ << ": ModelTrigger." << G4endl;

    bool istrigged = false;

    // only select the secondaries
    const G4Track* track = aFastTrack.GetPrimaryTrack();
    // secondaries
    if (track->GetParentID() != 0) {
        istrigged = true;
    }
    return istrigged;
}

void DummyFastSimG4Model::DoIt(const G4FastTrack& aFastTrack, G4FastStep& aFastStep) {
    // G4cout << __FILE__ << __LINE__ << ": DoIt." << G4endl;

    aFastStep.ProposeTrackStatus(fStopAndKill);
}
