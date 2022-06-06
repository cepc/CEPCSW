//////////////////////////////////////////////////////////////////
///
/// This is an interface of to call genfit
/// A genfit hit can be created
///
/// In this file, including:
///   a genfit hit class
///
///   Units are following GenfitUnit
#ifndef RECGENFITALG_GENFITHIT_H
#define RECGENFITALG_GENFITHIT_H
#include "TVector3.h"

namespace edm4hep{
    class SimTrackerHit;
    class TrackerHit;
}
namespace dd4hep {
    namespace DDSegmentation{
        class GridDriftChamber;
        class BitFieldCoder;
    }
}


class GenfitHit{
    public:
        GenfitHit(const edm4hep::TrackerHit* trackerHit,
                const edm4hep::SimTrackerHit* simTrackerHit,
                const dd4hep::DDSegmentation::BitFieldCoder* decoder,
                const dd4hep::DDSegmentation::GridDriftChamber* gridDriftChamber,
                double driftVelocity,double driftDistanceErr);
        ~GenfitHit(){;}
        unsigned long long getCellID()const;
        int getLayer()const;
        int getCell()const;
        double getDriftDistance()const{return m_driftDistance;}
        double getDriftDistanceErr()const{return m_driftDistanceErr;}
        double getDriftDistanceTruth()const{return m_driftDistanceTruth;}
        const edm4hep::SimTrackerHit* getSimTrackerHit()const{return m_simTrackerHit;}
        const edm4hep::TrackerHit* getTrackerHit()const{return m_trackerHit;}
        TVector3 getEnd0()const;
        TVector3 getEnd1()const;
        TVector3 getTruthPos()const;
        TVector3 getTruthMom()const;
        double getMaxDistance()const{return 0.6*1.4;}//FIXME
        int getLeftRightAmbig()const;

        void setDriftDistance(double d){m_driftDistance=d;}
        void setDriftDistanceErr(double de){m_driftDistanceErr=de;}
        void setDriftDistanceTruth(double dt){m_driftDistanceTruth=dt;}
        void print()const;

    private:
        const dd4hep::DDSegmentation::BitFieldCoder* m_decoder;
        const dd4hep::DDSegmentation::GridDriftChamber* m_gridDriftChamber;
        const edm4hep::TrackerHit* m_trackerHit;
        const edm4hep::SimTrackerHit* m_simTrackerHit;
        double m_driftDistance;
        double m_driftDistanceErr;
        double m_driftDistanceTruth;
        double m_driftVelocity;
};

#endif
