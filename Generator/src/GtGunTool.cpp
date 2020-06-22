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
    /*
    if (m_energies.value().size() != m_particles.value().size()) {
        error() << "Mismatched energies and particles." << endmsg;
        return StatusCode::FAILURE;
    }
    */
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

        TParticlePDG* particle = db_pdg->GetParticle(particle_name.c_str());
        if (particle) {
            pdgcode = particle->PdgCode();
            mass = particle->Mass(); // GeV
        } else {
            // guess it is pdg code
            pdgcode = atol(particle_name.c_str());
            if (!pdgcode) {
                error() << "Unsupported particle name/pdgcode " << particle_name << endmsg;
                return false;
            }
        }

        //double energy = m_energies.value()[i];
        double energy_min = m_energies_min.value()[0];
        double energy_max = m_energies_max.value()[0];

        // create the MC particle
        edm4hep::MCParticle mcp = event.m_mc_vec.create();
        mcp.setPDG(pdgcode);
        mcp.setGeneratorStatus(1);
        mcp.setSimulatorStatus(1);
        // mcp.setCharge();
        mcp.setTime(0.0);
        mcp.setMass(mass);
        // mcp.setVertex(); 
        // mcp.setEndpoint();

        // assume energy is momentum
        double p = energy_min==energy_max ? energy_max : CLHEP::RandFlat::shoot(energy_min, energy_max);
        double p1 = m_p1.value() ;
        
        // direction
        // by default, randomize the direction
        double theta = m_thetamins.value()[0]==m_thetamaxs.value()[0] ? m_thetamins.value()[0] : CLHEP::RandFlat::shoot(m_thetamins.value()[0], m_thetamaxs.value()[0]);
        double phi =   m_phimins  .value()[0]==m_phimaxs  .value()[0] ? m_phimins  .value()[0] : CLHEP::RandFlat::shoot(m_phimins  .value()[0], m_phimaxs  .value()[0]);
        double theta1 = m_dtheta.value() + theta;
        double phi1   = m_dphi.value()   + phi  ;
        double costheta = cos(theta*acos(-1)/180);
        double costheta1 = cos(theta1*acos(-1)/180);
        double phi_  = phi*acos(-1)/180;
        double phi1_ = phi1*acos(-1)/180;
        double sintheta = sqrt(1.-costheta*costheta);
        double sintheta1 = sqrt(1.-costheta1*costheta1);
        // check if theta min/max is set
        /*
        if (i < m_thetamins.value().size() 
            && i < m_thetamaxs.value().size()) {
            double thetamin = m_thetamins.value()[i];
            double thetamax = m_thetamaxs.value()[i];

            if (thetamin == thetamax) { // fixed theta
                costheta = cos(thetamin);
                sintheta = sin(thetamin);
                info() << "theta is fixed: " << thetamin << endmsg;
            }
        }

        if (i < m_phimins.value().size()
            && i < m_phimaxs.value().size()) {
            double phimin = m_phimins.value()[i];
            double phimax = m_phimaxs.value()[i];

            if (phimin == phimax) { // fixed phi
                phi = phimin;
                info() << "phi is fixed: " << phimin << endmsg;
            }
        }

        debug() << "Direction: "
                << " cos(theta): " << costheta
                << " phi: " << phi
                << endmsg;
        */ 
        double px = p*sintheta*cos(phi_);
        double py = p*sintheta*sin(phi_);
        double pz = p*costheta;
        std::cout<<"GenGt p="<<p<<", px="<<px<<",py="<<py<<",pz="<<pz<<",theta="<<theta<<",phi="<<phi<<std::endl;
        mcp.setMomentum(edm4hep::Vector3f(px,py,pz));
        // mcp.setMomentumAtEndpoint();
        // mcp.setSpin();
        // mcp.setColorFlow();
        if(m_create_second)
        {
            edm4hep::MCParticle mcp1 = event.m_mc_vec.create();
            mcp1.setPDG(pdgcode);
            mcp1.setGeneratorStatus(1);
            mcp1.setSimulatorStatus(1);
            mcp1.setTime(0.0);
            mcp1.setMass(mass);
            double px1 = p1*sintheta1*cos(phi1_);
            double py1 = p1*sintheta1*sin(phi1_);
            double pz1 = p1*costheta1;
            std::cout<<"GenGt p1="<<p1<<", px1="<<px1<<",py1="<<py1<<",pz1="<<pz1<<",theta1="<<theta1<<",phi1="<<phi1<<std::endl;
            mcp1.setMomentum(edm4hep::Vector3f(px1,py1,pz1));
        }

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

