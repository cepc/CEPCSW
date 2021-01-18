///////////////////////////////////////////////////////////
////   An implementation of genfit AbsBField - GenfitField
//// Authors:
////   Yao ZHANG (zhangyao@ihep.ac.cn)
///////////////////////////////////////////////////////////

#ifndef RECGENFITALG_GENFITFIELD_H
#define RECGENFITALG_GENFITFIELD_H

#include "AbsBField.h"
#include "DD4hep/Fields.h"

class TVector3;


/// Calss for translating Field into genfit::AbsBField
class GenfitField : public genfit::AbsBField{
    public:
        GenfitField(dd4hep::OverlayedField dd4hepField);
        virtual ~GenfitField(){;}

        //Get field value from TVector3
        TVector3 get(const TVector3& pos) const override;

        //Get field value from doubles
        void get(const double& posX, const double& posY, const double& posZ,
                double& Bx, double& By, double& Bz) const override;

        //Get Bz, kilogauss
        double getBz(const TVector3& pos) const;

        //Get dd4hep field
        const dd4hep::OverlayedField field() const {return m_dd4hepField;}

    private:
        dd4hep::OverlayedField m_dd4hepField;
};
#endif
