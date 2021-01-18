#include "BetheBlochEquationDedxSimTool.h"
#include "G4Step.hh"
#include "G4SystemOfUnits.hh"
#include "CLHEP/Random/RandGauss.h"

// https://folk.uib.no/ruv004/
DECLARE_COMPONENT(BetheBlochEquationDedxSimTool)


double BetheBlochEquationDedxSimTool::dedx(const G4Step* aStep)
{
    G4Track* gTrack = aStep->GetTrack() ;
    G4int z = gTrack->GetDefinition()->GetPDGCharge();
    if (z == 0) return 0;

    G4Material* material = gTrack->GetMaterial();
    G4double material_density = material->GetDensity() / (CLHEP::g/CLHEP::cm3); // conert from G4 unit.
    G4double material_Z = m_material_Z;
    G4double material_A = m_material_A;

    m_I = material_Z*10;  // Approximate

    G4double M = gTrack->GetDefinition()->GetPDGMass()/(CLHEP::eV);
    G4double gammabeta=aStep->GetPreStepPoint()->GetBeta() * aStep->GetPreStepPoint()->GetGamma();
    if(gammabeta<0.01)return 0;//too low momentum
    float beta = gammabeta/sqrt(1.0+pow(gammabeta,2));
    float gamma = gammabeta/beta;
    float Tmax = 2*m_me*pow(gammabeta,2)/(1+(2*gamma*m_me/M)+pow(m_me/M,2));
    float dedx = m_K*pow(z,2)*material_Z*(0.5*log(2*m_me*pow(gammabeta,2)*Tmax/pow(m_I,2))-pow(beta,2))/(material_A*pow(beta,2));    
    dedx = dedx*CLHEP::RandGauss::shoot(m_scale, m_resolution);
    double Dedx = dedx * (CLHEP::MeV/CLHEP::cm) ;
    return Dedx;
}

double BetheBlochEquationDedxSimTool::dedx(const edm4hep::MCParticle& mcp) {

    int z = mcp.getCharge();
    if (z == 0) return 0;
    float m_I = m_material_Z*10;  // Approximate

    double M = mcp.getMass();
    double gammabeta=sqrt(mcp.getMomentum()[0]*mcp.getMomentum()[0]+mcp.getMomentum()[1]*mcp.getMomentum()[1]+mcp.getMomentum()[2]*mcp.getMomentum()[2])/M;
    if(gammabeta<0.01)return 0;//too low momentum
    float beta = gammabeta/sqrt(1.0+pow(gammabeta,2));
    float gamma = gammabeta/beta;
    M = M*pow(10,9);//to eV
    float Tmax = 2*m_me*pow(gammabeta,2)/(1+(2*gamma*m_me/M)+pow(m_me/M,2));
    float dedx = m_K*pow(z,2)*m_material_Z*(0.5*log(2*m_me*pow(gammabeta,2)*Tmax/pow(m_I,2))-pow(beta,2))/(m_material_A*pow(beta,2));    
    dedx = dedx*CLHEP::RandGauss::shoot(m_scale, m_resolution);
    return dedx; // (CLHEP::MeV/CLHEP::cm) / (CLHEP::g/CLHEP::cm3)
}
double BetheBlochEquationDedxSimTool::dndx(double betagamma) {
    return -1;
}

StatusCode BetheBlochEquationDedxSimTool::initialize()
{
    m_me = 0.511*pow(10,6);//0.511 MeV to eV
    m_K = 0.307075;//const


    info() << "Initialize BetheBlochEquationDedxSimTool with following parameters" << endmsg;
    info() << "-> m_me: " << m_me << endmsg;
    info() << "-> m_K: " << m_K << endmsg;

    return StatusCode::SUCCESS;
}

StatusCode BetheBlochEquationDedxSimTool::finalize()
{
    return StatusCode::SUCCESS;
}
