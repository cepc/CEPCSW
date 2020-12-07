#include "DummyFastSimG4Model.h"

DummyFastSimG4Model::DummyFastSimG4Model(G4String aModelName, G4Region* aEnvelope)
    : G4VFastSimulationModel(aModelName, aEnvelope) {

}

DummyFastSimG4Model::~DummyFastSimG4Model() {

}

G4bool DummyFastSimG4Model::IsApplicable(const G4ParticleDefinition& aParticle) {
    return true;
}

G4bool DummyFastSimG4Model::ModelTrigger(const G4FastTrack& aFastTrack) {
    return true;
}

void DummyFastSimG4Model::DoIt(const G4FastTrack& aFastTrack, G4FastStep& aFastStep) {

}
