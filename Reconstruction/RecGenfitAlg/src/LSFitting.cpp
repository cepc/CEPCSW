#include "LSFitting.h"
#include <iostream>
#include <TMinuit.h>


std::vector<double> LSFitting::m_wireX;
std::vector<double> LSFitting::m_wireY;
std::vector<double> LSFitting::m_driftDist;

void LSFitting::Fitting(double& fitXc, double& fitYc, double& fitRadius)
{
    TMinuit* tMinuit = new TMinuit(3);  //initialize TMinuit with a maximum of 3 params
    tMinuit->SetFCN(FCN);
    double arglist[2];
    arglist[0] = 0;
    int ierflg = 0;
    tMinuit->SetPrintLevel(-1); // no print
    tMinuit->mnexcm("SET NOW",  arglist,1,ierflg); // no warning
    arglist[0] = 1;
    tMinuit->mnexcm("SET ERR", arglist ,1,ierflg);
    double xcErr= 1;
    double ycErr= 1;
    double radiusErr= 5;
    //fitXc=0;
    //fitYc=0;
    //fitRadius=1000;
    //std::cout<<" test Fitting---------------"<<fitXc<<" "<<fitYc<<" "<<fitRadius<<std::endl;
    tMinuit->mnparm(0,"xc",fitXc,xcErr/1.e4,fitXc-xcErr,fitXc+xcErr,ierflg);
    tMinuit->mnparm(1,"yc",fitYc,ycErr/1.e4,fitYc-ycErr,fitYc+ycErr,ierflg);
    tMinuit->mnparm(2,"r",fitRadius,radiusErr/1.e4,fitRadius-radiusErr,fitRadius+radiusErr,ierflg);
    arglist[0] = 0.0;
    arglist[1] = 1.0;
    tMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);
    double temp;
    tMinuit->GetParameter(0, fitXc,temp);
    tMinuit->GetParameter(1, fitYc,temp);
    tMinuit->GetParameter(2, fitRadius,temp);
    //std::cout<<" test Fitting---------------"<<fitXc<<" "<<fitYc<<" "<<fitRadius<<std::endl;

    delete tMinuit;
    return;
}

void LSFitting::FCN(int &npar, double *gin, double &f, double *par, int iflag)
{
    double fitXc=par[0];
    double fitYc=par[1];
    double fitRadius=par[2];
    f = 0.0;
    double hitChi2=0;

    for(unsigned int i=0;i<m_wireX.size();i++) {
        double fitDoca=fabs(sqrt((fitXc-m_wireX[i])*(fitXc-m_wireX[i])+
                    (fitYc-m_wireY[i])*(fitYc-m_wireY[i]))-fitRadius);
        double err=0.11;//mm
        hitChi2=(fitDoca-m_driftDist[i])*(fitDoca-m_driftDist[i])/(err*err);
        //std::cout<<"No."<<i<<" fitXY "<<fitXc<<" "<<fitYc<<" "<<fitRadius<<" wire "<<m_wireX[i]
        //    <<" "<<m_wireY[i]<<" hit, drift: "<<m_driftDist[i]<<", doca: "
            //<<fitDoca<<", hitchi: "<<hitChi2<<std::endl;
        f+= hitChi2;
        //std::cout<<"f: "<<f<<std::endl;
    }
    f = f/m_wireX.size();
    //std::cout<<"This fit has "<<m_wireX.size()<<" hits, chisq: "<<f<<std::endl;
}

void LSFitting::clear(){
    m_wireX.clear();
    m_wireY.clear();
    m_driftDist.clear();
    std::vector<double> tmp1;
    m_wireX.swap(tmp1);
    std::vector<double> tmp2;
    m_wireY.swap(tmp2);
    std::vector<double> tmp3;
    m_driftDist.swap(tmp3);
}

void LSFitting::setWire(double x,double y){
    m_wireX.push_back(x);
    m_wireY.push_back(y);
}
void LSFitting::setDrift(double driftDist){
    m_driftDist.push_back(driftDist);
}

void LSFitting::print(){
    std::cout<<" nHit "<<m_wireX.size()<<std::endl;
    for(unsigned int i=0;i<m_wireX.size();i++) {
        std::cout<<" wireX "<<m_wireX[i]<<" wireY "<<m_wireY[i]<<" drift "<<m_driftDist[i]<<std::endl;
    }
}
