#include "PrimaryGeneratorAction.h"

#include "G4Event.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4ParticleDefinition.hh"


PrimaryGeneratorAction::PrimaryGeneratorAction(ToolHandle<IG4PrimaryCnvTool>& cnvtool) 
    : G4VUserPrimaryGeneratorAction(), 
      tool(cnvtool) {

}

PrimaryGeneratorAction::~PrimaryGeneratorAction() {

}

void
PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent) {
    // Generate Vertex (G4PrimaryVertex) and Particle (G4PrimaryParticle).
    if (tool) {
        tool->mutate(anEvent);
    }

    // // Following is an example:
    // double x = 0.0;
    // double y = 0.0;
    // double z = 0.0;
    // double t = 0.0;
    // G4PrimaryVertex* g4vtx = new G4PrimaryVertex(x, y, z, t);


    // G4int pdgcode = 22;
    // // check the pdgid
    // G4ParticleTable* particletbl = G4ParticleTable::GetParticleTable();
    // G4ParticleDefinition* particle_def = particletbl->FindParticle(pdgcode);

    // double px = 0.0;
    // double py = 0.0;
    // double pz = 0.0;
    // G4PrimaryParticle* g4prim=new G4PrimaryParticle(particle_def, px, py, pz);
    // g4vtx->SetPrimary(g4prim);

    // anEvent->AddPrimaryVertex(g4vtx);
}

