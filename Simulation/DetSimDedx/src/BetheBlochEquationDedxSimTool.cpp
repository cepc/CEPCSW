#include "BetheBlochEquationDedxSimTool.h"
#include "G4Step.hh"

// https://folk.uib.no/ruv004/
DECLARE_COMPONENT(BetheBlochEquationDedxSimTool)


double BetheBlochEquationDedxSimTool::dedx(const G4Step* aStep)
{
    G4Track* gTrack = aStep->GetTrack() ;
    G4int z = gTrack->GetDefinition()->GetPDGCharge();
    if (z == 0) return 0;
    G4double M = gTrack->GetDefinition()->GetPDGMass();//MeV
    M = pow(10,6)*M; //to eV
    G4double gammabeta=aStep->GetPreStepPoint()->GetBeta() * aStep->GetPreStepPoint()->GetGamma();
    if(gammabeta<0.01)return 0;//too low momentum
    float beta = gammabeta/sqrt(1.0+pow(gammabeta,2));
    float gamma = gammabeta/beta;
    float Tmax = 2*m_me*pow(gammabeta,2)/(1+(2*gamma*m_me/M)+pow(m_me/M,2));
    float dedx = m_K*pow(z,2)*m_material_Z*(0.5*log(2*m_me*pow(gammabeta,2)*Tmax/pow(m_I,2))-pow(beta,2))/(m_material_A*pow(beta,2));    
    dedx = dedx*m_scale;// the material density can be absorbed in scale 
    dedx = dedx*(1+((*m_distribution)(m_generator)));
    return dedx*m_material_density; // MeV / cm 
}

StatusCode BetheBlochEquationDedxSimTool::initialize()
{
    m_distribution = new std::normal_distribution<double>(0, m_resolution);
    m_me = 0.511*pow(10,6);//0.511 MeV to eV
    m_K = 0.307075;//const
    m_I = m_material_Z*10;  // Approximate
    return StatusCode::SUCCESS;
}

StatusCode BetheBlochEquationDedxSimTool::finalize()
{
    return StatusCode::SUCCESS;
}
