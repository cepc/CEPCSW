#include "DataHelper/TrackHelper.h"
#include "DD4hep/DD4hepUnits.h"
#include "DataHelper/HelixClass.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include <iostream>

//get track position and momentum from TrackState
void CEPC::getPosMomFromTrackState(const edm4hep::TrackState& trackState,
        double Bz, TVector3& pos,TVector3& mom,double& charge,
        TMatrixDSym& covMatrix_6){
    double D0=trackState.D0;
    double phi=trackState.phi;
    double omega=trackState.omega;
    double Z0=trackState.Z0;
    double tanLambda=trackState.tanLambda;

    HelixClass helixClass;
    helixClass.Initialize_Canonical(phi,D0,Z0,omega,tanLambda,Bz);
    for(int i=0;i<3;i++){
        pos[i]=helixClass.getReferencePoint()[i];
        mom[i]=helixClass.getMomentum()[i];
    }

    charge = helixClass.getCharge();
    TMatrixDSym covMatrix_5(5);
    ///< lower triangular covariance matrix of the track parameters.
    ///  the order of parameters is  d0, phi, omega, z0, tan(lambda).
    std::array<float,15> covMatrix=trackState.covMatrix;
    int k=0;
    for(int i=0;i<5;i++){
        for(int j=0;j<5;j++){
            if(i>=j) { covMatrix_5(i,j)=covMatrix[k++]; }
        }
    }

    ///Error propagation
    //V(Y)=S * V(X) * ST , mS = S , mVy = V(Y) , helix.covariance() = V(X)
    TMatrix mS(covMatrix_6.GetNrows(),covMatrix_5.GetNrows());
    mS.Zero();
    double FCT = 2.99792458E-4;
    mS[0][0]=-sin(phi);
    mS[0][1]=-1*D0*cos(phi);
    mS[1][0]=cos(phi);
    mS[1][1]=-1*D0*sin(phi);
    mS[2][3]=1;
    mS[3][1]=FCT*Bz*(1/omega)*sin(phi);
    mS[3][2]=charge*FCT*Bz*(1/(omega*omega))*cos(phi);
    mS[4][1]=-1*FCT*Bz*(1/omega)*cos(phi);
    mS[4][2]=charge*FCT*Bz*(1/(omega*omega))*sin(phi);
    mS[5][2]=charge*tanLambda*Bz*(1/(omega*omega));
    mS[5][4]=-FCT*Bz/omega*(1+tanLambda*tanLambda);

    covMatrix_6= covMatrix_5.Similarity(mS);

    const char* m_name="TrackHelper";
    int m_debug=0;
    if(m_debug>=2){
        std::cout<<m_name<<" covMatrix_5 " <<std::endl;
        covMatrix_5.Print();
        std::cout<<m_name<<" mS " <<std::endl;
        mS.Print();
        std::cout<<m_name<<" covMatrix_6 " <<std::endl;
        covMatrix_6.Print();
        std::cout<<m_name<<" pos " <<std::endl;
        pos.Print();
        std::cout<<m_name<<" mom " <<std::endl;
        mom.Print();
    }
}

//Set track state from position, momentum and charge
void CEPC::getTrackStateFromPosMom(edm4hep::TrackState& trackState,double Bz,
        TVector3 pos,TVector3 mom,double charge,TMatrixDSym covMatrix_6){
    HelixClass helix;
    double pos_t[3]={pos.X(),pos.Y(),pos.Z()};
    double mom_t[3]={mom.X(),mom.Y(),mom.Z()};
    helix.Initialize_VP(pos_t,mom_t,charge,Bz);

    trackState.D0=helix.getD0();
    trackState.phi=helix.getPhi0();
    trackState.omega=helix.getOmega();
    trackState.Z0=helix.getZ0();
    trackState.tanLambda=helix.getTanLambda();
    trackState.referencePoint=helix.getReferencePoint();
    int m_debug=0;
    if(m_debug>=2){
        std::cout<<__FILE__<<" pos mom"<<std::endl;
        pos.Print();
        mom.Print();
        std::cout<<__FILE__<<" getTrackStateFromPosMom trackState  "<<trackState<<std::endl;
        std::cout<<__FILE__<<" getTrackStateFromPosMom cov6"<<std::endl;
        covMatrix_6.Print();
    }

    TMatrix Jacobian_matrix(5,6);
    Jacobian_matrix.Zero();

    double pt=sqrt(mom.X()*mom.X()+mom.Y()*mom.Y());
    Jacobian_matrix[0][0] = 2*pos.X()/(sqrt(pos.X()*pos.X()+pos.Y()*pos.Y()));
    Jacobian_matrix[0][1] = 2*pos.Y()/(sqrt(pos.X()*pos.X()+pos.Y()*pos.Y()));
    Jacobian_matrix[1][3] = mom.Y()/(pt*pt);
    Jacobian_matrix[1][4] = mom.X()/(pt*pt);
    Jacobian_matrix[2][3] = 2*charge*mom.X()/pt;
    Jacobian_matrix[2][4] = 2*charge*mom.Y()/pt;
    Jacobian_matrix[3][2] = 1.;
    Jacobian_matrix[4][3] = 2*mom.Z()*mom.X()/pt;
    Jacobian_matrix[4][4] = 2*mom.Z()*mom.Y()/pt;
    Jacobian_matrix[4][5] = 1./pt;

    TMatrixDSym covMatrix_5 = covMatrix_6.Similarity(Jacobian_matrix);

    std::array<float,15> covMatrix;
    int k=0;
    for(int i=0;i<5;i++){
        for(int j=0;j<5;j++){
            if(i>=j) { covMatrix[k++]=covMatrix_5(i,j); }
        }
    }
    trackState.covMatrix = covMatrix;
}
