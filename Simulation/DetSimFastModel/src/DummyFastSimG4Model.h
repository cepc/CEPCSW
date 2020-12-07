#ifndef DummyFastSimG4Model_h
#define DummyFastSimG4Model_h

#include "G4VFastSimulationModel.hh"

class DummyFastSimG4Model: public G4VFastSimulationModel {
public:

    DummyFastSimG4Model(G4String aModelName, G4Region* aEnvelope);
    ~DummyFastSimG4Model();

    virtual G4bool IsApplicable( const G4ParticleDefinition& aParticle );
    virtual G4bool ModelTrigger( const G4FastTrack& aFastTrack );
    virtual void DoIt( const G4FastTrack& aFastTrack, G4FastStep& aFastStep );

};

#endif

