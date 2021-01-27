//////////////////////////////////////////////////////////////////
///
/// This is an interface of to call genfit
/// A genfit track can be created, fitted and extrapolated
/// Track with only track representation(s)
///
/// In this file, including:
///   a genfit track class
///
///   Units are following DD4hepUnits
///
/// Authors:
///   Zhang Yao (zhangyao@ihep.ac.cn)
///   Y.Fujii (yfujii@ihep.ac.cn)
///   Yohei Nakatsugawa (yohei@ihep.ac.cn)
///
//////////////////////////////////////////////////////////////////

#ifndef RECGENFITALG_GENFITTRACK_H
#define RECGENFITALG_GENFITTRACK_H

#include "GenfitFitter.h"

//ROOT
#include "TVectorD.h"
#include "TMatrixDSym.h"

//STL
#include <vector>

class TLorentzVector;

namespace genfit{
    class Track;
    class FitStatus;
    class AbsTrackRep;
    class RKTrackRep;
    class KalmanFittedStateOnPlane;
}
namespace edm4hep{
    class MCParticle;
    class SimTrackerHitCollection;
    class ReconstructedParticle;
    class MCRecoTrackerAssociationCollection;
    class Track;
    class ConstTrackerHit;
    class Vector3d;
    class Vector3f;
}
namespace dd4hep {
    namespace DDSegmentation{
        class GridDriftChamber;
    }
}

class GenfitTrack {
    friend int GenfitFitter::processTrack(
            GenfitTrack* track, bool resort);

    friend int GenfitFitter::processTrackWithRep(
            GenfitTrack* track, int repID, bool resort);

    friend double GenfitFitter::extrapolateToHit(TVector3& poca,
            TVector3& pocaDir,
            TVector3& pocaOnWire, double& doca, const GenfitTrack* track,
            TVector3 pos, TVector3 mom, TVector3 end0, TVector3 end1, int debug,
            int repID, bool stopAtBoundary, bool calcJacobianNoise);

    friend double GenfitFitter::extrapolateToCylinder(TVector3& pos,
            TVector3& mom,
            GenfitTrack* track, double radius, const TVector3 linePoint,
            const TVector3 lineDirection, int hitID, int repID,
            bool stopAtBoundary, bool calcJacobianNoise);

    friend double GenfitFitter::extrapolateToPoint(TVector3& pos, TVector3& mom,
            const GenfitTrack* genfitTrack, const TVector3& point, int repID,
            bool stopAtBoundary, bool calcJacobianNoise) const;

    public:
    GenfitTrack(const GenfitField* field,
            const dd4hep::DDSegmentation::GridDriftChamber* seg,
            const char* name="GenfitTrack");
    virtual ~GenfitTrack();

    /// Add a Genfit track
    virtual bool createGenfitTrack(int pdgType,int charge,TLorentzVector pos, TVector3 mom,
            TMatrixDSym covM);
    //virtual bool createGenfitTrack(TLorentzVector posInit,TVector3 momInit,
            //TMatrixDSym covMInit);

    ///Create genfit track from MCParticle
    bool createGenfitTrackFromMCParticle(int pidTyep,const edm4hep::MCParticle&
            mcParticle, double eventStartTime=0.);
    bool createGenfitTrackFromEDM4HepTrack(int pidType,const edm4hep::Track& track,
            double eventStartTime);

    //  /// Prepare a hit list, return number of hits on track
    //  int PrepareHits();//TODO

    /// Add a space point measurement, return number of hits on track
    bool addSpacePointTrakerHit(edm4hep::ConstTrackerHit& hit, int hitID);

    /// Add a space point measurement, return number of hits on track
    virtual bool addSpacePointMeasurement(const TVectorD&, double,
            int detID=-1, int hitID=-1, bool smear=false);

    /// Add a WireMeasurement with MC truth position smeared by sigma
    virtual void addWireMeasurement(double driftDistance,
            double driftDistanceError, const TVector3& endPoint1,
            const TVector3& endPoint2, int lrAmbig, int detID, int hitID);

    /// Add a WireMeasurement with DC digi
    virtual bool addWireMeasurementOnTrack(edm4hep::Track& track, double sigma);

    ///Add space point from truth to track
    int addSimTrackerHits(const edm4hep::Track& track,
        const edm4hep::MCRecoTrackerAssociationCollection* assoHits,
        float sigma,bool smear=false);// float nSigmaSelection

    ///Store track to ReconstructedParticle
    bool storeTrack(edm4hep::ReconstructedParticle& dcRecParticle,int pidType,
            int ndfCut=1e9, double chi2Cut=1.e9);

    ///A tool to convert track to the first layer of DC
    void pivotToFirstLayer(edm4hep::Vector3d& pos,edm4hep::Vector3f& mom,
            edm4hep::Vector3d& firstPos, edm4hep::Vector3f& firstMom);

    /// Copy a track to event
    //void CopyATrack()const;

    ///Extrapolate to Hit
    double extrapolateToHit( TVector3& poca, TVector3& pocaDir,
            TVector3& pocaOnWire, double& doca, TVector3 pos, TVector3 mom,
            TVector3 end0,//one end of the hit wire
            TVector3 end1,//the orhter end of the hit wire
            int debug,
            int repID,
            bool stopAtBoundary,
            bool calcJacobianNoise) const;

    bool getPosMomCovMOP(int hitID, TLorentzVector& pos, TVector3& mom,
            TMatrixDSym& cov, int repID=0) const;

    /// get the seed position and momentum from track
    const TLorentzVector getSeedStatePos() const;
    const TVector3 getSeedStateMom() const;
    void getSeedStateMom(TLorentzVector& pos, TVector3& mom) const;
    unsigned int getNumPoints() const;

    /// get the fitted track status
    const genfit::FitStatus* getFitStatus(int repID=0) const;
    int getFittedState(TLorentzVector& pos, TVector3& mom, TMatrixDSym& cov,
            int repID=0, bool biased=true) const;
    int getNumPointsWithFittedInfo(int repID=0) const;
    bool getFirstPointWithFittedInfo(int repID=0) const;
    bool fitSuccess(int repID) const;

    /// get the wire infomation
    int getDetIDWithFitterInfo(int hitID, int idRaw=0) const;

    int getPDG(int id=0) const;
    int getPDGCharge(int id=0) const;

    /// print genfit track
    void printSeed() const;//print seed
    void printFitted(int repID=0) const;//print seed
    void print(TLorentzVector pos, TVector3 mom, const char* comment="") const;

    /// set and get debug level
    void setDebug(int debug);
    int getDebug(void) const { return m_debug;}

    /// get name of this object
    const char* getName() const {return m_name;}
    genfit::Track* getTrack() const{return m_track;}

    /// Add a track representation
    genfit::RKTrackRep* addTrackRep(int pdg);

    protected:
    //genfit::Track* getTrack() {return m_track;}

    private:

    const char* m_name;

    ///Note! private functions are using genfit unit, cm and MeV

    genfit::AbsTrackRep* getRep(int id=0) const {return m_reps[id];}
    bool getMOP(int hitID, genfit::MeasuredStateOnPlane& mop,
            genfit::AbsTrackRep* trackRep=nullptr) const;

    genfit::Track* m_track;/// track
    std::vector<genfit::AbsTrackRep*> m_reps;/// track repesentations
    int m_debug;/// debug level

    const GenfitField* m_genfitField;//pointer to genfit field
    const dd4hep::DDSegmentation::GridDriftChamber* m_gridDriftChamber;

    static const int s_PDG[2][5];

};

#endif
