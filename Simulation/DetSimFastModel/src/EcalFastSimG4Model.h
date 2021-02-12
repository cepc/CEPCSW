#ifndef EcalFastSimG4Model_h
#define EcalFastSimG4Model_h

#include "G4VFastSimulationModel.hh"

class EcalFastSimG4Model: public G4VFastSimulationModel {
public:

    EcalFastSimG4Model(G4String aModelName, G4Region* aEnvelope);
    ~EcalFastSimG4Model();

    virtual G4bool IsApplicable( const G4ParticleDefinition& aParticle );
    virtual G4bool ModelTrigger( const G4FastTrack& aFastTrack );
    virtual void DoIt( const G4FastTrack& aFastTrack, G4FastStep& aFastStep );

};

#endif

