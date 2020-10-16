#include "G4PrimaryCnvTool.h"

#include "G4Event.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4ParticleDefinition.hh"

DECLARE_COMPONENT(G4PrimaryCnvTool)

bool G4PrimaryCnvTool::mutate(G4Event* anEvent) {

    auto mcCol = m_mcParCol.get();
    info() << "Start a new event: " << endmsg;
    for ( auto p : *mcCol ) {
        info() << p.getObjectID().index << " : [";
        for ( auto it = p.daughters_begin(), end = p.daughters_end(); it != end; ++it ) {
            info() << " " << it->getObjectID().index;
        }
        info() << " ]; " << endmsg;

        // only the GeneratorStatus == 1 is used.
        if (p.getGeneratorStatus() != 1) {
            continue;
        }

        // vertex
        const edm4hep::Vector3d& vertex = p.getVertex();
        double t = p.getTime()*CLHEP::ns;
        G4PrimaryVertex* g4vtx = new G4PrimaryVertex(vertex.x*CLHEP::mm,
                                                     vertex.y*CLHEP::mm,
                                                     vertex.z*CLHEP::mm,
                                                     t);

        info() << "Geant4 Primary Vertex: ("
               << vertex.x*CLHEP::mm << ","
               << vertex.y*CLHEP::mm << ","
               << vertex.z*CLHEP::mm << ")"
               << endmsg;

        // pdg/particle
        int pdgcode = p.getPDG();
        G4ParticleTable* particletbl = G4ParticleTable::GetParticleTable();
        G4ParticleDefinition* particle_def = particletbl->FindParticle(pdgcode);

        // momentum
        const edm4hep::Vector3f& momentum = p.getMomentum();
        G4PrimaryParticle* g4prim = new G4PrimaryParticle(particle_def,
                                                          momentum.x*CLHEP::GeV,
                                                          momentum.y*CLHEP::GeV,
                                                          momentum.z*CLHEP::GeV);

        g4vtx->SetPrimary(g4prim);

        anEvent->AddPrimaryVertex(g4vtx);
    }

    return true;
}
