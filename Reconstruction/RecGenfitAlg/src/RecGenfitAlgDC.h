//////////////////////////////////////////////////////////////////////
///
/// This is an algorithm for track fitting for CEPC track with genfit.
///
/// In this file, including:
///   An algorithm for track fitting with genfit with 5 pid hypothesis
///
///   Units are following DD4hepUnits
///
/// Authors:
///   Yao ZHANG(zhangyao@ihep.ac.cn)
///
/////////////////////////////////////////////////////////////////////

#ifndef RECGENFITALG_RECGENFITALGDC_H
#define RECGENFITALG_RECGENFITALGDC_H

#include "GaudiAlg/GaudiAlgorithm.h"
#include "GaudiKernel/NTuple.h"
#include "k4FWCore/DataHandle.h"
#include "DD4hep/Fields.h"
#include <string>

class GenfitFitter;
class GenfitField;
class GenfitTrack;
class IGeomSvc;
class time;
namespace genfit{
    class EventDisplay;
}
namespace dd4hep {
    class Detector;
    //class rec::CellIDPositionConverter;
    namespace DDSegmentation{
        class GridDriftChamber;
        class BitFieldCoder;
    }
}
namespace edm4hep{
    class EventHeaderCollection;
    class MCParticleCollection;
    class SimTrackerHitCollection;
    class TrackCollection;
    class TrackerHitCollection;
    class MCRecoTrackerAssociationCollection;
    class ReconstructedParticle;
    class ReconstructedParticleCollection;
}

/////////////////////////////////////////////////////////////////////////////

class RecGenfitAlgDC:public GaudiAlgorithm {
    public:
        RecGenfitAlgDC (const std::string& name, ISvcLocator* pSvcLocator);
        StatusCode initialize() override;
        StatusCode execute() override;
        StatusCode finalize() override;

    private:
        GenfitFitter* m_genfitFitter;//The pointer to a GenfitFitter
        const GenfitField* m_genfitField;//The pointer to a GenfitField

        StatusCode fitMC();
        void debugTrack(int pidType,const GenfitTrack* genfitTrack);
        void debugEvent();
        //void debugDoca(RecMdcTrack* recMdcTrack,const GenfitTrack* genfitTrack);
        //void debugTruthHit();

        DataHandle<edm4hep::EventHeaderCollection> _headerCol{
            "EventHeaderCol", Gaudi::DataHandle::Reader, this};
        //Drift chamber rec hit and trac
        DataHandle<edm4hep::TrackerHitCollection> m_digiDCHitsCol{
            "DigiDCHitsCollection", Gaudi::DataHandle::Reader, this};
        DataHandle<edm4hep::TrackCollection> m_dcTrackCol{
            "DCTrackCollection", Gaudi::DataHandle::Reader, this};
        //Mc truth
        DataHandle<edm4hep::MCParticleCollection> m_mcParticleCol{
            "MCParticle", Gaudi::DataHandle::Reader, this};
        DataHandle<edm4hep::SimTrackerHitCollection> m_simDCHitCol{
            "DriftChamberHitsCollection" , Gaudi::DataHandle::Reader, this};
        //output hits and particles
        DataHandle<edm4hep::TrackerHitCollection> m_dcFitRecHitCol{
            "DCFitRecHitsCollection", Gaudi::DataHandle::Writer, this};
        DataHandle<edm4hep::ReconstructedParticleCollection> m_dcRecParticleCol{
            "DCRecParticleCollection", Gaudi::DataHandle::Writer, this};

        const unsigned int m_nPDG;//5:e,mu,pi,K,proton
        SmartIF<IGeomSvc> m_geomSvc;
        dd4hep::OverlayedField m_dd4hepField;
        dd4hep::Detector* m_dd4hep;
        dd4hep::DDSegmentation::GridDriftChamber* m_gridDriftChamber;
        dd4hep::DDSegmentation::BitFieldCoder* m_decoder;
        Gaudi::Property<std::string> m_readout_name{this, "readout",
            "DriftChamberHitsCollection"};
        int m_inputType;
        int m_debug;
        std::string m_fitterType;//default DAFRef, DAFRef, DAF
        bool m_correctBremsstrahlung;//defalut do not correct bremsstrahlung
        bool m_noMaterialEffects;//defalut false, correct material effects
        int m_maxIteration;//default 20
        bool m_resortHits;//default true
        double m_bStart;
        double m_bFinal;
        double m_dcCornerCuts;
        int m_ndfCut;//default true
        double m_chi2Cut;//default 1000
        int m_fitSuccess[5];
        int m_nDCTrack;
        int m_debugPid;
        bool m_useTruthSeed;
        bool m_useRecLRAmbig;

        std::string m_genfitHistRootName;
        bool m_firstPid;
        bool m_useTruthHit;
        bool m_showDisplay;
        double m_initCovResPos;
        double m_initCovResMom;
        genfit::EventDisplay* m_genfitDisplay;
        clock_t m_timer;

        /// tuples
        NTuple::Tuple*  m_tuple;
        NTuple::Item<int> m_run;
        NTuple::Item<int> m_evt;
        NTuple::Item<int> m_tkId;
        NTuple::Item<int> m_mcIndex;//number of navigated mcParicle
        NTuple::Matrix<double> m_pocaPosMc;//2 dim matched particle and 3 pos.
        NTuple::Matrix<double> m_pocaMomMc;//2 dim matched particle and 3 mom.
        NTuple::Array<double> m_pocaMomMcP;//2 dim matched particle and p
        NTuple::Array<double> m_pocaMomMcPt;//2 dim matched particle and pt
        NTuple::Array<double> m_pocaPosMdc;//pos 0:x,1:y,2:z
        NTuple::Array<double> m_pocaMomMdc;//mom. 0:px,1:py,2:pz
        NTuple::Item<int> m_pidIndex;
        NTuple::Matrix<double> m_firstPosKal;//5 hyposis and pos. at first
        NTuple::Array<double> m_firstMomKalP;//5 hyposis and mom. at first
        NTuple::Array<double> m_firstMomKalPt;//5 hyposis and mom. at first
        NTuple::Matrix<double> m_pocaPosKal;//5 hyposis and 3 mom.
        NTuple::Matrix<double> m_pocaMomKal;//5 hyposis and 3 mom.
        NTuple::Array<double> m_pocaMomKalP;//5 hyposis and p
        NTuple::Array<double> m_pocaMomKalPt;//5 hyposis and pt
        NTuple::Array<int> m_nDofKal;
        NTuple::Array<double> m_chi2Kal;
        NTuple::Array<bool> m_convergeKal;
        NTuple::Array<bool> m_isFittedKal;
        NTuple::Item<int> m_nDigi;
        NTuple::Item<int> m_nHitMc;
        NTuple::Item<int> m_nHitMdc;
        NTuple::Item<int> m_nClusterCgem;
        NTuple::Item<int> m_nHitKalInput;
        NTuple::Array<double> m_mdcHitDriftT;
        NTuple::Array<double> m_mdcHitDriftDl;
        NTuple::Array<double> m_mdcHitDriftDr;
        NTuple::Array<int> m_mdcHitLr;
        NTuple::Array<int> m_mdcHitLayer;
        NTuple::Array<int> m_mdcHitWire;
        NTuple::Array<double> m_mdcHitExpDoca;
        NTuple::Array<double> m_mdcHitExpMcDoca;
        NTuple::Array<double> m_mdcHitErr;
        NTuple::Array<int> m_nHitFailedKal;
        NTuple::Array<int> m_nHitFitted;
        NTuple::Array<double> m_time;
        //truth
        NTuple::Array<int> m_mdcHitMcLr;
        NTuple::Array<int> m_mdcHitMcTkId;
        NTuple::Array<double> m_mdcHitMcDrift;
        NTuple::Array<double> m_mdcHitMcX;
        NTuple::Array<double> m_mdcHitMcY;
        NTuple::Array<double> m_mdcHitMcZ;
        NTuple::Array<double> m_mdcHitExpMcPocaX;
        NTuple::Array<double> m_mdcHitExpMcPocaY;
        NTuple::Array<double> m_mdcHitExpMcPocaZ;
        NTuple::Array<double> m_mdcHitExpMcPocaWireX;
        NTuple::Array<double> m_mdcHitExpMcPocaWireY;
        NTuple::Array<double> m_mdcHitExpMcPocaWireZ;

};
#endif
