#include "GtGunTool.h"

#include "TDatabasePDG.h"
#include "CLHEP/Units/SystemOfUnits.h"
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGauss.h"

DECLARE_COMPONENT(GtGunTool)

StatusCode
GtGunTool::initialize() {
    StatusCode sc;
    // the particles names/pdgs and energies should be specified.
    if (m_particles.value().size()==0) {
        error() << "Please specify the list of particle names/pdgs" << endmsg;
        return StatusCode::FAILURE;
    }

    // Position
    if (m_positionXs.value().size()
        && m_positionXs.value().size() != m_particles.value().size()) {
        error() << "Mismatched PositionXs and particles." << endmsg;
        return StatusCode::FAILURE;
    }
    if (m_positionYs.value().size()
        && m_positionYs.value().size() != m_particles.value().size()) {
        error() << "Mismatched PositionYs and particles." << endmsg;
        return StatusCode::FAILURE;
    }
    if (m_positionZs.value().size()
        && m_positionZs.value().size() != m_particles.value().size()) {
        error() << "Mismatched PositionZs and particles." << endmsg;
        return StatusCode::FAILURE;
    }
    
    // Energy
    if (m_energymins.value().size() != m_particles.value().size()) {
        error() << "Mismatched energies and particles." << endmsg;
        return StatusCode::FAILURE;
    }
    
    // others should be empty or specify
    if (m_thetamins.value().size()
        && m_thetamins.value().size() != m_particles.value().size()) {
        error() << "Mismatched thetamins and particles." << endmsg;
        return StatusCode::FAILURE;
    }
    if (m_thetamaxs.value().size()
        && m_thetamaxs.value().size() != m_particles.value().size()) {
        error() << "Mismatched thetamaxs and particles." << endmsg;
        return StatusCode::FAILURE;
    }

    if (m_phimins.value().size()
        && m_phimins.value().size() != m_particles.value().size()) {
        error() << "Mismatched phimins and particles." << endmsg;
        return StatusCode::FAILURE;
    }
    if (m_phimaxs.value().size()
        && m_phimaxs.value().size() != m_particles.value().size()) {
        error() << "Mismatched phimaxs and particles." << endmsg;
        return StatusCode::FAILURE;
    }

    return sc;
}

StatusCode
GtGunTool::finalize() {
    StatusCode sc;
    return sc;
}

bool
GtGunTool::mutate(MyHepMC::GenEvent& event) {

    TDatabasePDG* db_pdg = TDatabasePDG::Instance();

    // The default unit here is GeV.
    // but we don't add the unit, because in geant4 it is multiplied.

    for (int i = 0; i < m_particles.value().size(); ++i) {
        const std::string& particle_name = m_particles.value()[i];
        int pdgcode = 0;
        double mass = 0;
        double charge = 0;

        TParticlePDG* particle = db_pdg->GetParticle(particle_name.c_str());
        if (particle) {
            pdgcode = particle->PdgCode();
            mass = particle->Mass(); // GeV
            charge = particle->Charge()/3; // in ROOT, it is in units of |e|/3
        } else {
            // guess it is pdg code
            pdgcode = atol(particle_name.c_str());
            if (!pdgcode) {
                error() << "Unsupported particle name/pdgcode " << particle_name << endmsg;
                return false;
            }
        }

        double energy = m_energymins.value()[i]==m_energymaxs.value()[i] ? m_energymins.value()[i] : CLHEP::RandFlat::shoot(m_energymins.value()[i], m_energymaxs.value()[i]);

        // create the MC particle
        edm4hep::MCParticle mcp = event.m_mc_vec.create();
        mcp.setPDG(pdgcode);
        mcp.setGeneratorStatus(1);
        mcp.setSimulatorStatus(1);
        mcp.setCharge(static_cast<float>(charge));
        mcp.setTime(0.0);
        mcp.setMass(mass);

        // Unit is mm
        double x = 0;
        double y = 0;
        double z = 0;
        if (i<m_positionXs.value().size()) { x = m_positionXs.value()[i]; }
        if (i<m_positionYs.value().size()) { y = m_positionYs.value()[i]; }
        if (i<m_positionZs.value().size()) { z = m_positionZs.value()[i]; }

        mcp.setVertex(edm4hep::Vector3d(x,y,z)); 
        // mcp.setEndpoint();

        // assume energy is momentum
        double p = energy;
        
        // direction
        // by default, randomize the direction
        double theta = m_thetamins.value()[i]==m_thetamaxs.value()[i] ? m_thetamins.value()[i] : CLHEP::RandFlat::shoot(m_thetamins.value()[i], m_thetamaxs.value()[i]);
        double phi =   m_phimins  .value()[i]==m_phimaxs  .value()[i] ? m_phimins  .value()[i] : CLHEP::RandFlat::shoot(m_phimins  .value()[i], m_phimaxs  .value()[i]);
        double costheta = cos(theta*acos(-1)/180);
        double phi_  = phi*acos(-1)/180;
        double sintheta = sqrt(1.-costheta*costheta);
        double px = p*sintheta*cos(phi_);
        double py = p*sintheta*sin(phi_);
        double pz = p*costheta;
        std::cout<<"GenGt p="<<p<<", px="<<px<<",py="<<py<<",pz="<<pz<<",theta="<<theta<<",phi="<<phi<<std::endl;
        mcp.setMomentum(edm4hep::Vector3f(px,py,pz));
        // mcp.setMomentumAtEndpoint();
        // mcp.setSpin();
        // mcp.setColorFlow();

    }

    // event.SetEventHeader( m_processed_event, -99, 9999, "Generator");

    return true;
}

bool
GtGunTool::finish() {
    return true;
}

bool
GtGunTool::configure_gentool() {

    return true;
}

