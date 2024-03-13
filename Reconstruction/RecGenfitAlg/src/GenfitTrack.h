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
///
//////////////////////////////////////////////////////////////////

#ifndef RECGENFITALG_GENFITTRACK_H
#define RECGENFITALG_GENFITTRACK_H

#include "GenfitFitter.h"
#include "GenfitHit.h"

//ROOT
#include "TVector3.h"
#include "TVectorD.h"
#include "TMatrixDSym.h"
#include "TLorentzVector.h"

//Gaudi
#include "GaudiKernel/SmartIF.h"

//STL
#include <vector>

class TLorentzVector;
class IGeomSvc;
class WireMeasurementDC;

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
    class MutableReconstructedParticle;
    class MCRecoTrackerAssociationCollection;
    class Track;
    class MutableTrack;
    class TrackCollection;
    class TrackerHit;
    class SimTrackerHit;
    class Vector3d;
    class Vector3f;
    class TrackerHitCollection;
}
namespace dd4hep {
    namespace DDSegmentation{
        class GridDriftChamber;
        class BitFieldCoder;
    }
    namespace rec{
        class ISurface;
    }
}

class GenfitTrack {
    friend int GenfitFitter::processTrack(
            GenfitTrack* track, bool resort);
    
    friend int GenfitFitter::processTrackWithRep(
            GenfitTrack* track, int repID, bool resort);

    public:
    GenfitTrack(const GenfitField* field,
            const dd4hep::DDSegmentation::GridDriftChamber* seg,
            SmartIF<IGeomSvc> geom,
            const char* name="GenfitTrack");
    virtual ~GenfitTrack();

    ///Create genfit track from MCParticle
    bool createGenfitTrackFromMCParticle(int pidTyep,const edm4hep::MCParticle&
            mcParticle, double eventStartTime=0.);
    bool createGenfitTrackFromEDM4HepTrack(int pidType,const edm4hep::Track& track,
            double eventStartTime,bool isUseCovTrack);

    /// ---------Add measurements---------
    ///Add one space point measurement, return number of hits on track
    virtual bool addSpacePointMeasurement(const TVector3&,std::vector<float>
        sigmaU,std::vector<float> sigmaV,int cellID,int hitID);

    ///Add silicon space points from edm4hep::track
    int addSpacePointsSi(const edm4hep::Track& track,
            std::vector<float> sigmaU,std::vector<float> sigmaV);

    ///Add drift chamber space points from edm4hep::track
    int addSpacePointsDC(const edm4hep::Track& track,
            const edm4hep::MCRecoTrackerAssociationCollection* assoHits,
        std::vector<float> sigmaU,std::vector<float> sigmaV);

    ///Add WireMeasurements of hits on track
    virtual int addWireMeasurementsOnTrack(edm4hep::Track& track,float sigma,
            const edm4hep::MCRecoTrackerAssociationCollection* assoHits,
            int sortMethod,bool truthAmbig,float skipCorner, float skipNear);

    ///Add WireMeasurements of hits on track from hit selection
    virtual int addWireMeasurementsFromList(std::vector<edm4hep::TrackerHit*>& hits,float sigma,
            const edm4hep::MCRecoTrackerAssociationCollection* assoHits,
            int sortMethod,bool truthAmbig,float skipCorner, float skipNear);

    virtual int addWireMeasurementsFromListTrF(const edm4hep::TrackerHitCollection* trkHits,
            float sigma,int sortMethod);

    ///Add one silicon hits
    bool addSiliconMeasurement(edm4hep::TrackerHit* hit,
            float sigmaU,float sigmaV,int cellID,int hitID);

    ///Add silicon measurements, return number of hits on track
    int addSiliconMeasurements(edm4hep::Track& track,
            std::vector<float> sigmaU,std::vector<float> sigmaV);

    bool debugDistance(const edm4hep::TrackerHitCollection* dCDigiCol,
            int& nFittedDC,int& nFittedSDT,int& ngenfitHit,
            std::vector<double>& smearDistance,
            std::vector<double>& truthDistance,double driftVelocity);
    bool GetDocaRes(int hitID,  double& DriftDis,double& fittedDoca,
            double& res,int repID=0, bool biased=true) const;

    ///Store track to ReconstructedParticle
    bool storeTrack(edm4hep::MutableReconstructedParticle& recParticle,
            edm4hep::MutableTrack& track,
            TVector3& pocaToOrigin_pos,
            TVector3& pocaToOrigin_mom,
            TMatrixDSym& pocaToOrigin_cov,
            int pidType, int ndfCut, double chi2Cut,
            int& nFittedDC, int& nFittedSDT,int& ngenfitHit,
            std::vector<double>& trackL, std::vector<double>& hitMom,
            std::vector<float>& truthMomEdep,
            const edm4hep::MCRecoTrackerAssociationCollection* assoHits,
            std::vector<double>& driftDis,
            std::vector<double>& FittedDoca,
            std::vector<double>& truthDoca,
            std::vector<double>& Res,
            std::vector<double>& truthRes);

    ///A tool to convert track to the first layer of DC
  
    void pivotToFirstLayer(const edm4hep::Vector3d& pos,
            const edm4hep::Vector3f& mom, edm4hep::Vector3d& firstPos,
            edm4hep::Vector3f& firstMom);

    /// Copy a track to event
    //void CopyATrack()const;

    //return dd4hep unit
    double extrapolateToHit(TVector3& poca, TVector3& pocaDir,
            TVector3& pocaOnWire, double& doca,edm4hep::MCParticle mcParticle,
            int cellID, int repID, bool stopAtBoundary, bool calcJacobianNoise)const;
    /// Extrapolate the track to the point
    /// Output: pos and mom of POCA point to point
    /// Input: genfitTrack,point,repID,stopAtBoundary and calcAverageState
    /// repID same with pidType
    double extrapolateToPoint(TVector3& pos, TVector3& mom, TMatrixDSym& cov,
            const TVector3& point, int repID=0,
            bool stopAtBoundary = false, bool calcJacobianNoise = true) const;

    /// Extrapolate the track to the cyliner at fixed raidus
    /// Output: pos and mom at fixed radius
    /// Input: genfitTrack, radius of cylinder at center of the origin,
    ///        repID, stopAtBoundary and calcAverageState
    double extrapolateToCylinder(TVector3& pos, TVector3& mom,
            double radius, const TVector3 linePoint,
            const TVector3 lineDirection, int hitID =0, int repID=0,
            bool stopAtBoundary=false, bool calcJacobianNoise=true);


    bool getPosMomCovMOP(int hitID, TLorentzVector& pos, TVector3& mom,
            TMatrixDSym& cov, int repID=0) const;

    /// get the seed position and momentum from track
    const TLorentzVector getSeedStatePos() const;
    const TVector3 getSeedStateMom() const;
    void getSeedStateMom(TLorentzVector& pos, TVector3& mom) const;
    unsigned int getNumPoints() const;
    unsigned int getNumPointsDet(int cellID) const;

    /// get the fitted track status
    const genfit::FitStatus* getFitStatus(int repID=0) const;
    int getFittedState(TLorentzVector& pos, TVector3& mom, TMatrixDSym& cov,
            int trackPointId=0, int repID=0, bool biased=true) const;
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
    void setDebug(int debug){m_debug=debug;}
    void setDebugGenfit(int debug);
    void setDebugLocal(int debug);
    int getDebug(void) const { return m_debug;}

    /// get name of this object
    const char* getName() const {return m_name;}
    genfit::Track* getTrack() const{return m_track;}

    /// Add a track representation
    genfit::RKTrackRep* addTrackRep(int pdgType,int charge);

    /// Get a hit according to index
    GenfitHit* GetHit(long unsigned int i) const;
    
    static void getTrackFromEDMTrackFinding(const edm4hep::Track& edm4HepTrack,
            double& charge, TVectorD& trackParam, TMatrixDSym& cov,TVector3& pos,
            TVector3& mom);
    protected:
    //genfit::Track* getTrack() {return m_track;}

    private:

    /// ---------Add a Genfit track-------
    bool createGenfitTrack(int pdgType,int charge,
            TVectorD trackParam, TMatrixDSym covMInit_6);

    int getDetTypeID(unsigned long long cellID) const;
    const char* m_name;

    ///Note! private functions are using genfit unit, cm and MeV

    genfit::AbsTrackRep* getRep(int id=0) const;
    bool getMOP(int hitID, genfit::MeasuredStateOnPlane& mop,
            genfit::AbsTrackRep* trackRep=nullptr) const;
    const dd4hep::rec::ISurface* getISurface(edm4hep::TrackerHit* hit);
    void getSeedCov(TMatrixDSym& cov);
    void getAssoSimTrackerHit(
            const edm4hep::MCRecoTrackerAssociationCollection*& assoHits,
            edm4hep::TrackerHit* trackerHit,
            edm4hep::SimTrackerHit& simTrackerHit) const;
    void getEndPointsOfWire(int cellID,TVector3& end0,TVector3& end1)const;
    void getTrackFromEDMTrack(const edm4hep::Track& edm4HepTrack,
            double& charge, TVectorD& trackParam, TMatrixDSym& cov) const;
    void getTrackFromMCPartile(const edm4hep::MCParticle mcParticle,
            TVectorD& trackParam, TMatrixDSym& cov) const;
    void getPosMomFromMCPartile(const edm4hep::MCParticle mcParticle,
            TVector3& pos,TVector3& mom) const;
    void clearGenfitHitVec();
    void getISurfaceOUV(const dd4hep::rec::ISurface* iSurface,TVector3& o,
            TVector3& u,TVector3& v,double& lengthU,double& lengthV);
    void getMeasurementAndCov(edm4hep::TrackerHit* hit,TVector3& pos,TMatrixDSym& cov);
    int getSigmas(int cellID,std::vector<float> sigmaUVec,
        std::vector<float> sigmaVVec,float& sigmaU,float& sigmaV)const;
    bool isCDCHit(edm4hep::TrackerHit* hit);
    GenfitHit* makeAGenfitHit(edm4hep::TrackerHit* trackerHit,
            edm4hep::SimTrackerHit* simTrackerHitAsso,
            double sigma,bool truthAmbig,double skipCorner,double skipNear);
    void getSortedTrackerHits(std::vector<edm4hep::TrackerHit*>& hits,
            const edm4hep::MCRecoTrackerAssociationCollection* assoHits,
            std::vector<edm4hep::TrackerHit*>& sortedDCTrackerHits,
            int sortMethod);
    void getSortedTrackerHitsTrF(std::vector<edm4hep::TrackerHit*> trackerHits,
            std::vector<edm4hep::TrackerHit*>& sortedDCTrackerHits,
            int sortMethod=1);

    genfit::Track* m_track;/// track
    int m_debug;/// debug level
    int m_debugLocal;/// debug level local

    SmartIF<IGeomSvc> m_geomSvc;
    const GenfitField* m_genfitField;//pointer to genfit field
    const dd4hep::DDSegmentation::GridDriftChamber* m_gridDriftChamber;
    const dd4hep::DDSegmentation::BitFieldCoder* m_decoderDC;

    static const int s_PDG[2][5];

    std::vector<GenfitHit*> m_genfitHitVec;

};
#endif
