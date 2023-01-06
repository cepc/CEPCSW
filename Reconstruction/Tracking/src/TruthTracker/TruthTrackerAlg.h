#ifndef TruthTrackerAlg_h
#define TruthTrackerAlg_h

#include "GaudiAlg/GaudiAlgorithm.h"
#include "k4FWCore/DataHandle.h"
#include "DD4hep/Fields.h"
#include "GaudiKernel/NTuple.h"

#include "TRandom3.h"

class IGeomSvc;
namespace dd4hep {
    class Detector;
    namespace DDSegmentation{
        class GridDriftChamber;
        class BitFieldCoder;
    }
}
namespace edm4hep {
    class MCParticleCollection;
    class SimTrackerHitCollection;
    class TrackerHitCollection;
    class TrackCollection;
    class Track;
    class MutableTrack;
    class TrackState;
    class ReconstructedParticleCollection;
    class MCRecoTrackerAssociationCollection;
    class MCRecoParticleAssociationCollection;
}

class TruthTrackerAlg: public GaudiAlgorithm
{
    public:
        TruthTrackerAlg(const std::string& name, ISvcLocator* svcLoc);

        virtual StatusCode initialize() override;
        virtual StatusCode execute() override;
        virtual StatusCode finalize() override;

    private:
        void getTrackStateFromMcParticle(const edm4hep::MCParticleCollection*
                mcParticleCol, edm4hep::TrackState& stat);
        int addSimHitsToTk(DataHandle<edm4hep::SimTrackerHitCollection>&
                colHandle, edm4hep::TrackerHitCollection*& truthTrackerHitCol,
                edm4hep::MutableTrack& track, const char* msg,int nHitAdded);
        int smearDCTkhit(DataHandle<edm4hep::TrackerHitCollection>&
                colHandle,DataHandle<edm4hep::TrackerHitCollection>& smearCol,
                DataHandle<edm4hep::SimTrackerHitCollection>& SimDCHitCol,
                DataHandle<edm4hep::SimTrackerHitCollection>& SimSmearDCHitCol,
                DataHandle<edm4hep::MCRecoTrackerAssociationCollection>& AssoDCHitCol,
                DataHandle<edm4hep::MCRecoTrackerAssociationCollection>& AssoSmearDCHitCol,
                double resX, double resY, double resZ);
        int addHitsToTk(DataHandle<edm4hep::TrackerHitCollection>&
                colHandle, edm4hep::MutableTrack& track, const char* msg,int nHitAdded);
        int addIdealHitsToTk(DataHandle<edm4hep::TrackerHitCollection>&
                colHandle, edm4hep::TrackerHitCollection*& truthTrackerHitCol,
                edm4hep::MutableTrack& track, const char* msg,int nHitAdded);

        int addHotsToTk(edm4hep::Track& sourceTrack,edm4hep::MutableTrack&
                targetTrack, int hitType,const char* msg,int nHitAdded);
        int nHotsOnTrack(edm4hep::Track& track, int hitType);
        int trackerHitColSize(DataHandle<edm4hep::TrackerHitCollection>& hitCol);
        int simTrackerHitColSize(DataHandle<edm4hep::SimTrackerHitCollection>& hitCol);
        bool getTrackStateFirstHit(
                DataHandle<edm4hep::SimTrackerHitCollection>& dcSimTrackerHitCol,
                float charge,edm4hep::TrackState& trackState);
        SmartIF<IGeomSvc> m_geomSvc;
        dd4hep::Detector* m_dd4hep;
        dd4hep::OverlayedField m_dd4hepField;
        dd4hep::DDSegmentation::GridDriftChamber* m_gridDriftChamber;
        dd4hep::DDSegmentation::BitFieldCoder* m_decoder;
        void debugEvent();
        //unit length is mm
        void getCircleFromPosMom(double pos[3],double mom[3],
                double Bz,double q,double& helixRadius,double& helixXC,double& helixYC);
        int makeNoiseHit(edm4hep::SimTrackerHitCollection* SimVec,
                edm4hep::TrackerHitCollection* Vec,
                edm4hep::MCRecoTrackerAssociationCollection* AssoVec,
                const edm4hep::TrackerHitCollection* digiDCHitsCol,
                const edm4hep::MCRecoTrackerAssociationCollection* assoHits);
        bool debugNoiseHitsCol(edm4hep::TrackerHitCollection* Vec);

            //reader
            //        DataHandle<edm4hep::TrackerHitCollection> m_NoiseHitCol{
            //            "NoiseDCHitsCollection", Gaudi::DataHandle::Reader, this};
        DataHandle<edm4hep::MCParticleCollection> m_mcParticleCol{
            "MCParticle", Gaudi::DataHandle::Reader, this};
        DataHandle<edm4hep::SimTrackerHitCollection> m_DCSimTrackerHitCol{
            "DriftChamberHitsCollection", Gaudi::DataHandle::Reader, this};
        DataHandle<edm4hep::TrackerHitCollection> m_DCDigiCol{
            "DigiDCHitCollection", Gaudi::DataHandle::Reader, this};
        DataHandle<edm4hep::MCRecoTrackerAssociationCollection>
            m_DCHitAssociationCol{ "DCHitAssociationCollection",
                Gaudi::DataHandle::Reader, this};
        DataHandle<edm4hep::TrackCollection>
            m_siSubsetTrackCol{ "SubsetTracks",
                Gaudi::DataHandle::Reader, this};
        DataHandle<edm4hep::TrackerHitCollection> m_SITSpacePointCol{
            "SITSpacePoints" , Gaudi::DataHandle::Reader, this};
        //        DataHandle<edm4hep::TrackerHitCollection> m_SETSpacePointCol{
        //            "SETSpacePoints" , Gaudi::DataHandle::Reader, this};
        DataHandle<edm4hep::TrackerHitCollection> m_FTDSpacePointCol{
            "FTDSpacePoints" , Gaudi::DataHandle::Reader, this};
        DataHandle<edm4hep::TrackerHitCollection> m_VXDTrackerHits{
            "VXDTrackerHits" , Gaudi::DataHandle::Reader, this};
        DataHandle<edm4hep::TrackerHitCollection> m_SETTrackerHits{
            "SETTrackerHits" , Gaudi::DataHandle::Reader, this};
        DataHandle<edm4hep::TrackerHitCollection> m_SITTrackerHits{
            "SITTrackerHits" , Gaudi::DataHandle::Reader, this};
        DataHandle<edm4hep::TrackerHitCollection> m_FTDTrackerHits{
            "FTDTrackerHits" , Gaudi::DataHandle::Reader, this};
        DataHandle<edm4hep::SimTrackerHitCollection> m_VXDCollection{
            "VXDCollection" , Gaudi::DataHandle::Reader, this};
        DataHandle<edm4hep::SimTrackerHitCollection> m_SETCollection{
            "SETCollection" , Gaudi::DataHandle::Reader, this};
        DataHandle<edm4hep::SimTrackerHitCollection> m_SITCollection{
            "SITCollection" , Gaudi::DataHandle::Reader, this};
        DataHandle<edm4hep::SimTrackerHitCollection> m_FTDCollection{
            "FTDCollection" , Gaudi::DataHandle::Reader, this};
        //writer
        DataHandle<edm4hep::TrackCollection> m_DCTrackCol{
            "DCTrackCollection", Gaudi::DataHandle::Writer, this};
        DataHandle<edm4hep::TrackCollection> m_SDTTrackCol{
            "SDTTrackCollection", Gaudi::DataHandle::Writer, this};
        DataHandle<edm4hep::TrackerHitCollection> m_truthTrackerHitCol{
            "TruthTrackerHitCollection", Gaudi::DataHandle::Writer, this};

        // Smear hit
        DataHandle<edm4hep::SimTrackerHitCollection> w_SimSmearHCol{
            "SmearSimHitsCollection", Gaudi::DataHandle::Writer, this};
        DataHandle<edm4hep::TrackerHitCollection> m_SmeartruthTrackerHitCol{
            "SmearTrackerHitCollection", Gaudi::DataHandle::Writer, this};
        DataHandle<edm4hep::MCRecoTrackerAssociationCollection> w_SmearAssociationCol{
            "SmearDCHitAssociationCollection", Gaudi::DataHandle::Writer, this};
        // make noise hit
        DataHandle<edm4hep::SimTrackerHitCollection> w_SimNoiseHCol{
            "NoiseSimHitsCollection", Gaudi::DataHandle::Writer, this};
        DataHandle<edm4hep::TrackerHitCollection> w_NoiseHitCol{
            "NoiseDCHitsCollection", Gaudi::DataHandle::Writer, this};
        DataHandle<edm4hep::MCRecoTrackerAssociationCollection> w_NoiseAssociationCol{
            "NoiseDCHitAssociationCollection", Gaudi::DataHandle::Writer, this};


        //readout for getting segmentation
        Gaudi::Property<std::string> m_readout_name{this, "readout",
            "DriftChamberHitsCollection"};
        Gaudi::Property<bool> m_hist{this,"hist",false};
        Gaudi::Property<bool> m_useDC{this,"useDC",true};
        Gaudi::Property<bool> m_useSi{this,"useSi",true};
        Gaudi::Property<bool> m_useTruthTrack{this,"useTruthTrack",false};
        Gaudi::Property<bool> m_useSiTruthHit{this,"useSiTruthHit",false};
        Gaudi::Property<bool> m_skipSecondaryHit{this,"skipSecondaryHit",true};
        Gaudi::Property<bool> m_useFirstHitForDC{this,"useFirstHitForDC",false};
        Gaudi::Property<bool> m_useSiSpacePoint{this,"useSiSpacePoint",false};
        Gaudi::Property<bool> m_useIdealHit{this,"useIdealHit",false};

        Gaudi::Property<bool> m_useNoiseHits{this,"useNoiseHits",false};
        Gaudi::Property<bool> m_smearHits{this,"smearHits",false};

        Gaudi::Property<float> m_momentumCut{this,"momentumCut",0.1};//momentum cut for the first hit
        Gaudi::Property<float> m_resPT{this,"resPT",0};//ratio
        Gaudi::Property<float> m_resPz{this,"resPz",0};//ratio
        Gaudi::Property<float> m_resX{this,"resX",0.11};//mm
        Gaudi::Property<float> m_resY{this,"resY",0.11};//mm
        Gaudi::Property<float> m_resZ{this,"resZ",0.11};//mm
        Gaudi::Property<float> m_driftVelocity{this,"driftVelocity",40};// um/us
        Gaudi::Property<float> m_resMomPhi{this,"resMomPhi",0};//radian
        Gaudi::Property<float> m_resMomTheta{this,"resMomTheta",0};//radian
        Gaudi::Property<float> m_resVertexX{this,"resVertexX",0.003};//3um
        Gaudi::Property<float> m_resVertexY{this,"resVertexY",0.003};//3um
        Gaudi::Property<float> m_resVertexZ{this,"resVertexZ",0.003};//3um
        Gaudi::Property<int> m_maxDCDigiCut{this,"maxDigiCut",1e6};
        Gaudi::Property<std::vector<float> > m_resVXD{this,"resVXD",{0.003,0.003,0.003}};//mm
        Gaudi::Property<std::vector<float> > m_resSIT{this,"resSIT",{0.003,0.003,0.003}};//mm
        Gaudi::Property<std::vector<float> > m_resSET{this,"resSET",{0.003,0.003,0.003}};//mm
        Gaudi::Property<std::vector<float> > m_resFTDPixel{this,"resFTDPixel",{0.003,0.003,0.003}};//mm
        Gaudi::Property<std::vector<float> > m_resFTDStrip{this,"resFTDStrip",{0.003,0.003,0.003}};//mm

        Gaudi::Property<double> m_fHitPurity{this,"fHitPurity",0.1};
        Gaudi::Property<float> m_pocaTime  { this, "pocaTime", 225};// ns

        double m_helixRadius,m_helixXC,m_helixYC;
        double m_helixRadiusFirst,m_helixXCFirst,m_helixYCFirst;

        NTuple::Tuple*  m_tuple;
        NTuple::Item<int> m_run;
        NTuple::Item<int> m_evt;
        NTuple::Array<double> m_siMom;
        NTuple::Array<double> m_siPos;
        NTuple::Array<double> m_mcMom;
        NTuple::Array<double> m_mcPos;

        NTuple::Item<int> m_nDCTrackHit;
        NTuple::Item<int> m_nSmearDCTrackHit;
        NTuple::Item<int> m_nNoiseDCTrackHit;
        NTuple::Array<float> m_DriftDistance;
        NTuple::Array<float> m_SmearDriftDistance;
        NTuple::Array<float> m_NoiseDriftDistance;

        NTuple::Item<int> m_nSimTrackerHitVXD;
        NTuple::Item<int> m_nSimTrackerHitSIT;
        NTuple::Item<int> m_nSimTrackerHitSET;
        NTuple::Item<int> m_nSimTrackerHitFTD;
        NTuple::Item<int> m_nSimTrackerHitDC;
        NTuple::Item<int> m_nTrackerHitVXD;
        NTuple::Item<int> m_nTrackerHitSIT;
        NTuple::Item<int> m_nTrackerHitSET;
        NTuple::Item<int> m_nTrackerHitFTD;
        NTuple::Item<int> m_nTrackerHitDC;
        NTuple::Item<int> m_nTrackerHitErrVXD;
        NTuple::Item<int> m_nTrackerHitErrSIT;
        NTuple::Item<int> m_nTrackerHitErrSET;
        NTuple::Item<int> m_nTrackerHitErrFTD;
        NTuple::Item<int> m_nSpacePointSIT;
        NTuple::Item<int> m_nSpacePointSET;
        NTuple::Item<int> m_nSpacePointFTD;
        NTuple::Item<int> m_nSpacePointErrVXD;
        NTuple::Item<int> m_nSpacePointErrSIT;
        NTuple::Item<int> m_nSpacePointErrSET;
        NTuple::Item<int> m_nSpacePointErrFTD;
        NTuple::Item<int> m_nHitOnSiTkVXD;
        NTuple::Item<int> m_nHitOnSiTkSIT;
        NTuple::Item<int> m_nHitOnSiTkSET;
        NTuple::Item<int> m_nHitOnSiTkFTD;
        NTuple::Item<int> m_nHitOnSdtTkVXD;
        NTuple::Item<int> m_nHitOnSdtTkSIT;
        NTuple::Item<int> m_nHitOnSdtTkSET;
        NTuple::Item<int> m_nHitOnSdtTkFTD;
        NTuple::Item<int> m_nHitOnSdtTkDC;
        NTuple::Item<int> m_nHitOnSdtTk;
        NTuple::Item<int> m_nNoiseOnSdtTk;

        TRandom3 fRandom;
};

#endif
