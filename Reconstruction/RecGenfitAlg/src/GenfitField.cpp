#include "GenfitField.h"
#include "GenfitUnit.h"

//External
#include "DD4hep/DD4hepUnits.h"

//genfit
#include "FieldManager.h"
#include "ConstField.h"

///-------------------------------------------------------------
/// GenfitField constructor
///-------------------------------------------------------------
GenfitField::GenfitField(dd4hep::OverlayedField dd4hepField)
:m_dd4hepField(dd4hepField)
{
    genfit::FieldManager::getInstance()->init(this);
}


///-------------------------------------------------------------
/// GenfitField get field value by TVector3
///-------------------------------------------------------------
/// unit in genfit is cm and kGauss
TVector3 GenfitField::get(const TVector3& pos) const
{
    double B[3]={1e9,1e9,1e9};
    get(pos.X(),pos.Y(),pos.Z(),B[0],B[1],B[2]);
    return TVector3(B[0],B[1],B[2]);
}

///-------------------------------------------------------------
/// GenfitField get field value by double
///-------------------------------------------------------------
/// unit in genfit for position is cm and Bxyz is kGauss
/// unit in DD4hepUnits kilogaus and cm? FIXME
void
GenfitField::get(const double& posX, const double& posY, const double& posZ,
        double& Bx, double& By, double& Bz) const
{
    /// get field from dd4hepField
    const dd4hep::Direction& field=m_dd4hepField.magneticField(
            dd4hep::Position(
              posX/GenfitUnit::cm*dd4hep::cm,
              posY/GenfitUnit::cm*dd4hep::cm,
              posZ/GenfitUnit::cm*dd4hep::cm));

    Bx=field.X()/dd4hep::kilogauss*GenfitUnit::kilogauss;
    By=field.Y()/dd4hep::kilogauss*GenfitUnit::kilogauss;
    Bz=field.Z()/dd4hep::kilogauss*GenfitUnit::kilogauss;

    //std::cout<<"GenfitField "
    //  <<Form("xyz(%f,%f,%f)cm B(%f,%f,%f)kilogauss",
    //      posX/GenfitUnit::cm*dd4hep::cm,
    //      posY/GenfitUnit::cm*dd4hep::cm,
    //      posZ/GenfitUnit::cm*dd4hep::cm,Bx,By,Bz)
    //  <<std::endl;
}

double GenfitField::getBz(const TVector3& pos) const
{
    return get(pos).Z();
}
