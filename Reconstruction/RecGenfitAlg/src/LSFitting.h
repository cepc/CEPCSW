#ifndef RECONSTRUCTION_RECGENFITALG_LSFITTING_H
#define RECONSTRUCTION_RECGENFITALG_LSFITTING_H
#include <vector>
#include <cmath>

class LSFitting{
public:
    void Fitting(double& fitXc, double& fitYc, double& fitRadisu);
    static void FCN(int &npar, double *gin, double &f, double *par, int iflag);
    void setWire(double x,double y);
    void setDrift(double driftDist);
    void print();
    void clear();
    static std::vector<double> m_wireX;
    static std::vector<double> m_wireY;
    static std::vector<double> m_driftDist;
};
#endif
