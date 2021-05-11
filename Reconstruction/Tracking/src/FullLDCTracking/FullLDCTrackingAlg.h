#ifndef FULLLDCTRACKING_H
#define FULLLDCTRACKING_H 1

#include "k4FWCore/DataHandle.h"
#include "GaudiAlg/GaudiAlgorithm.h"
#include "GearSvc/IGearSvc.h"
#include "TrackSystemSvc/ITrackSystemSvc.h"

#include "DataHelper/ClusterExtended.h"
#include "DataHelper/TrackExtended.h"
#include "DataHelper/TrackerHitExtended.h"
#include "DataHelper/TrackHitPair.h"
#include "DataHelper/HelixClass.h"
// #include "DataHelper/ClusterShapes.h"
#include "DataHelper/GroupTracks.h"

#include <map>
#include <set>

#include "TrackSystemSvc/MarlinTrkUtils.h"
#include "TrackSystemSvc/IMarlinTrack.h"

#include "edm4hep/SimTrackerHitCollection.h"
#include "edm4hep/TrackerHitCollection.h"
#include "edm4hep/TrackCollection.h"

#include <UTIL/BitField64.h>
#include <UTIL/ILDConf.h>

//using namespace edm4hep ;

//class gear::GearMgr ;

namespace MarlinTrk {
  class HelixFit;
  class IMarlinTrkSystem ;
}


/** === FullLDCTracking_MarlinTrk Processor === <br>
 * Processor performing track finding procedure in 
 * the entire LDC detector by linking track segments
 * found by the SiliconTracking module in the silicon detectors
 * and by the LEPTracking module in TPC. 
 * <h4>Input collections and prerequisites</h4> 
 * Processor requires collection of digitized vertex, sit, ftd, set, etd & tpc tracker hits 
 * and also the collections of tracks found in the silicon detectors
 * and in TPC.
 * <h4>Output</h4>
 * Processor produces an LCIO collection of the Tracks. Each track is characterised by 
 * five parameters : Omega (signed curvuture), Tan(lambda) where
 * lambda is the dip angle, Phi (azimuthal angle @ point of closest approach), D0 (signed impact parameter),
 * Z0 (displacement along z axis at the point of closest approach to IP). 
 * Covariance matrix for these parameters is also provided.
 * Only lower left corner of the covariance matrix is stored. The sequence of the covariance matrix elements 
 * assigned to track is the following: <br>
 * (D0,D0) <br>
 * (Phi,D0), (Phi,Phi) <br>
 * (Omega,D0), (Omega,Phi), (Omega,Omega) <br>
 * (Z0,D0), (Z0,Phi), (Z0,Omega), (Z0,Z0) <br>
 * (TanL,D0), (TanL,Phi), (TanL,Omega), (TanL,Z0), (TanL,TanL) <br>
 * The number of hits in the different subdetectors associated
 * with each track can be accessed via method Track::getSubdetectorHitNumbers().
 * This method returns vector of integers : <br>
 * number of VTX hits used in the track fit is the 1st element in this vector  
 * (Track::getSubdetectorHitNumbers()[0]) <br>
 * number of FTD hits used in the track fit is the 2nd element in this vector  
 * (Track::getSubdetectorHitNumbers()[1]) <br>
 * number of SIT hits used in the track fit is the 3d element in this vector  
 * (Track::getSubdetectorHitNumbers()[2]) <br>
 * number of TPC hits used in the track fit is the 4th element in this vector  
 * (Track::getSubdetectorHitNumbers()[3]) <br>
 * number of SET hits used in the track fit is the 5th element in this vector  
 * (Track::getSubdetectorHitNumbers()[4]) <br>
 * number of ETD hits used in the track fit is the 6th element in this vector  
 * (Track::getSubdetectorHitNumbers()[5]) <br>
 * total number of VTX hits in track is the 7th element in this vector 
 * (Track::getSubdetectorHitNumbers()[6]) <br>
 * total number of FTD hits in track is the 8th element in this vector
 * (Track::getSubdetectorHitNumbers()[7]) <br>
 * total number of SIT hits in track is the 9th element in this vector
 * (Track::getSubdetectorHitNumbers()[8]) <br>
 * total number of TPC hits in track is the 10th element in this vector
 * (Track::getSubdetectorHitNumbers()[9]) <br>
 * total number of SET hits in track is the 11th element in this vector
 * (Track::getSubdetectorHitNumbers()[10]) <br>
 * total number of ETD hits in track is the 12th element in this vector
 * (Track::getSubdetectorHitNumbers()[11]) <br>
 * Output track collection has by default a name "LDCTracks". 
 * @param VTXHitCollection name of input VTX TrackerHit collection <br>
 * (default parameter value : "VTXTrackerHits") <br>
 * @param FTDPixelHitCollectionName name of input FTD Pixel TrackerHit collection <br>
 * (default parameter value : "FTDPixelTrackerHits") <br>
 * @param FTDSpacePointCollectionName name of input FTD Space Point TrackerHit collection <br>
 * (default parameter value : "FTDSpacePoints") <br>
 * @param SITHitCollection name of input SIT TrackerHit collection <br>
 * (default parameter value : "SITTrackerHits") <br>
 * @param TPCHitCollection name of input TPC TrackerHit collection <br>
 * (default parameter value : "TPCTrackerHits") <br>
 * @param SETHitCollection name of input SET TrackerHit collection <br>
 * (default parameter value : "SETTrackerHits") <br>
 * @param ETDHitCollection name of input ETD TrackerHit collection <br>
 * (default parameter value : "ETDTrackerHits") <br>
 * @param TPCTracks collection name of TPC tracks <br>
 * (default parameter value : "TPCTracks") <br>
 * @param TPCTracksMCPRelColl Name of input TPC track to MC particle relation collection <br>
 * (default parameter value : "TPCTracksMCP") <br>
 * @param SiTracks collection name of Si tracks <br>
 * (default parameter value : "SiTracks") <br>
 * @param SiTracksMCPRelColl Name of input Si track to MC particle relation collection <br>
 * (default parameter value : "SiTracksMCP") <br> 
 * @param LDCTrackCollection name of the output LDC track collection <br>
 * (default parameter value : "LDCTracks") <br>
 * @param Chi2FitCut cut on the Chi2/Ndf of the track fit <br>
 * (default parameter value : 100.0) <br>
 * @param Chi2PrefitCut cut on the prefit Chi2 of the track candidate, 
 * prefit is done with the simple helix hypothesis <br>
 * (default parameter value : 1e+5) <br>
  * @param AngleCutForMerging  cut on opening angle between 
 * particle momentum reconstructed with TPC and momentum reconstructed
 * with the Silicon detectors. If the opening angle is smaller that this cut
 * the track segment in Silicon trackers and in TPC are tested for their
 * compatibility <br>
 * (default parameter value : 0.10) <br>
 * @param OmegaCutForMerging  cut on the relative difference in the track Omega
 * parameter reconstructed with TPC and with Si detectors. If the relative difference is smaller
 * than this cut, the track segments in TPC and Si are tested for their compatibility <br>
 * (default parameter value : 0.25) <br>
 * @param D0CutForMerging Upper cutoff on the difference in D0 [mm] to allow for merging 
 * of the Si and TPC segments <br>
 * (default parameter value : 500) <br>
 * @param Z0CutForMerging Upper cutoff on the difference in Z0 [mm] to allow for merging
 * of the Si and TPC segments <br>
 * (default parameter value : 1000) <br>
 * @param Debug flag to allow for printout of debug information,
 * if set to 1 debugging printout is activated
 * (default parameter value : 1) <br>
 * @param ForceSiTPCMerging This flag steers merging of Si and TPC track segments. If ForceMerging=1
 * Si and TPC track segments are forced to be merged if the opening angle between Si track 
 * momentum and TPC track momentum
 * is less than AngleCutForForcedMerging (see below) and difference in tracks 
 * parameters Omega is less than OmegaCutForForcedMerging (see below) <br>
 * (default parameter value : 0)
 * @param AngleCutForForcedMerging cut on opening angle between Si track momentum and
 * TPC track momentum. Used to steer forced merging of Si and TPC track segments. <br>
 * (default parameter value : 0.05)
 * @param OmegaCutForForcedMerging cut on the difference between Si and TPC tracks parameter
 * Omega. Used to steer forced merging of Si and TPC track segments. Relative 
 * errors are compared. <br>
 * (default parameter value : 0.15) <br>
 * @param D0CutForForcedMerging Upper cutoff on the difference in D0 to allow for forced
 * merging of the Si and TPC segments <br>
 * (default parameter value : 50) <br>
 * @param Z0CutForForcedMerging Upper cutoff on the difference in Z0 to allow for forced
 * merging of the Si and TPC segments <br>
 * (default parameter value : 200) <br>
 * @param ForceTPCSegmentsMerging If this flag is set to 1, the code attempts to 
 * merge TPC segments from the low pt splitted loopers <br>
 * (default parameter value : 1) <br>
 * @param D0CutToMergeTPCSegments cut on the difference in the track parameter
 * d0 [mm] to allow for merging TPC segments <br>
 * (default parameter value : 100) <br>
 * @param Z0CutToMergeTPCSegments cut on the difference in the track parameter
 * z0 [mm] to allow for merging TPC segments <br>
 * (default parameter value : 5000) <br> 
 * @param DeltaPCutToMergeTPCSegments cut on the magnitude [GeV/c] of the vectorial difference
 * of the momentum vectors, associated with TPC segments, for the TPC segment's merging procedure <br>
 * (default parameter value : 0.1) <br>
 * @param PtCutToMergeTPCSegments lower cutoff on Pt of the TPC segments of the looping track for
 * the merging procedure.
 * If transverse momentum of the segments is less than cutoff the segments are allowed to be merged. <br>
 * (default parameter value : 1.2) <br> 
 * @param AssignTPCHits If this flag is set to 1, the code attempts to assign left-over 
 * TPC hits to the accepted track candidates. No track refit is done in case when hit is assigned
 * to the existing track <br>
 * (default parameter value : 1) <br>
 * @param AssignETDHits If this flag is set to 1, the code attempts to assign  
 * ETD hits to the accepted track candidates. No track refit is done in case when hit is assigned
 * to the existing track <br>
 * (default parameter value : 1) <br>
 * @param AssignVTXHits If this flag is set to 1, the code attempts to assign left-over 
 * VTX hits to the accepted track candidates. Track refit is done in case when hit is assigned
 * to the existing track <br>
 * (default parameter value : 1) <br>
 * @param AssignFTDHits If this flag is set to 1, the code attempts to assign left-over 
 * FTD hits to the accepted track candidates. Track refit is done in case when hit is assigned
 * to the existing track <br>
 * (default parameter value : 1) <br>
 * @param AssignSITHits If this flag is set to 1, the code attempts to assign left-over 
 * SIT hits to the accepted track candidates. Track refit is done in case when hit is assigned
 * to the existing track <br>
 * (default parameter value : 1) <br>
 * @param AssignSETHits If this flag is set to 1, the code attempts to assign  
 * SET hits to the accepted track candidates. Track refit is done in case when hit is assigned
 * to the existing track <br>
 * (default parameter value : 1) <br>
 * @param TPCHitToTrackDistance Cut on the distance between left-over TPC hit and the track helix
 * to allow for assignment of the hit with a given track <br>
 * (default parameter value : 15.0) <br>
 * @param VTXHitToTrackDistance Cut on the distance between left-over VTX hit and the track helix
 * to allow for assignment of the hit with a given track <br>
 * (default parameter value : 1.5) <br>
 * @param FTDHitToTrackDistance Cut on the distance between left-over FTD hit and the track helix
 * to allow for assignment of the hit with a given track <br>
 * (default parameter value : 2.0) <br>
 * @param SITHitToTrackDistance Cut on the distance between left-over SIT hit and the track helix
 * to allow for assignment of the hit with a given track <br>
 * (default parameter value : 2.0) <br>
 * @param SETHitToTrackDistance Cut on the distance between SET hit and the track helix
 * to allow for assignment of the hit with a given track <br>
 * (default parameter value : 2.0) <br>
 * @param ETDHitToTrackDistance Cut on the distance between ETD hit and the track helix
 * to allow for assignment of the hit with a given track <br>
 * (default parameter value : 10.0) <br>
 * @param NHitsExtrapolation Number of the last track hits for extrapolating helix
 * to the outer tracking detectors (SET, ETD) <br>
 * (default parameter value : 35) <br>
 * @param CutOnTPCHits minimal number of TPC hits, used in the track fit, which is 
 * required for tracks which have no hits from the Si detectors <br>
 * (default parameter value : 35) <br> 
 * @param CutOnTrackD0 cut on the d0 parameter of the track. If the d0 parameter is greater that 
 * this cut, track is rejected <br>
 * (default parameter value : 500) <br>
 * @param CutOnTrackZ0 cut on the z0 parameter of the track. If the z0 parameter is greater that 
 * this cut, track is rejected <br>
 * (default parameter value : 500) <br>
 * @param ForbidOverlapInZTPC If this flag is set to 1 then merging of the TPC semiloops is 
 * forbiden for segment overlapping in z <br>
 * (default parameter value : 0) <br>
 * @param ForbidOverlapInZComb If this flag is set to 1 then merging of left-over TPC semiloop and
 * combined Si-TPC track is their segments overlap in z <br>
 * (default parameter value : 0) <br>
 * @param cosThetaCutHighPtMerge cut on cos theta between the two momentum vectors 
 * when considering merger of high Pt tracks <br>
 * (default is 0.99) <br>
 * @param cosThetaCutSoftHighPtMerge cut on the cos theta between the two momentum vectors 
 * when considering merger of high Pt tracks for softer dp/p cut <br>
 * (default is 0.998) <br>
 * @param momDiffCutHighPtMerge cut on dp/p 
 * when considering merger of high Pt tracks <br>
 * (default is 0.01 [1/GeV]) <br>
 * @param momDiffCutSoftHighPtMerge softer cut on dp/p  
 * when considering merger of high Pt tracks <br>
 * (default is 0.25 [1/GeV]) <br>
 * @param hitDistanceCutHighPtMerge cut on 3D distance between hit 
 * and helix extrapolation when considering merger of high Pt tracks <br>
 * (default is 25.0 [mm]) <br>
 * @param maxHitDistanceCutHighPtMerge cut for max 3D distance between any hit 
 * and helix extrapolation when considering merger of high Pt tracks <br>
 * (default is 50.0 [mm]) <br>
 * @param maxFractionOfOutliersCutHighPtMerge cut on maximum fraction of outliers 
 * when considering merger of high Pt tracks <br>
 * (default is 0.95 ) <br>
 
 
 * @author A. Raspereza (MPI Munich)<br>
 */

class FullLDCTrackingAlg : public GaudiAlgorithm {
  
public:
  
  FullLDCTrackingAlg(const std::string& name, ISvcLocator* svcLoc);
  virtual StatusCode initialize();
  virtual StatusCode execute();
  virtual StatusCode finalize();
  
protected:
  
  void prepareVectors();
  void CleanUp();
  void MergeTPCandSiTracks();
  void MergeTPCandSiTracksII();

  TrackExtended* CombineTracks(TrackExtended* tpcTrk, TrackExtended* siTrk, float maxAllowedOutliers, bool testCombinationOnly);

  // TrackExtended* TrialCombineTracks(TrackExtended* tpcTrk, TrackExtended* siTrk);

  void Sorting(TrackExtendedVec& trackVec);
  void SelectCombinedTracks();
  void AddNotCombinedTracks();
  void CheckTracks();
  void AddNotAssignedHits();
  // void RemoveSplitTracks();
  void AddTrackColToEvt(TrackExtendedVec & trkVec, edm4hep::TrackCollection* trkCol);
  float CompareTrk(TrackExtended * first, TrackExtended * second, float d0Cut, float z0Cut, int iopt);
  
  float CompareTrkII(TrackExtended * first, TrackExtended * second, float d0Cut, float z0Cut, int iopt, float &Angle);
  float CompareTrkIII(TrackExtended * first, TrackExtended * second, float d0Cut, float z0Cut, int iopt, float &Angle);
  
  void SortingTrackHitPairs(TrackHitPairVec & trackHitPairVec); 
  void AssignSiHitsToTracks(TrackerHitExtendedVec hitVec, float dcut);
  
  void AssignTPCHitsToTracks(TrackerHitExtendedVec hitVec, float dcut);
  
  void AssignOuterHitsToTracks(TrackerHitExtendedVec hitVec, float dcut, int refit);
  
  void CreateExtrapolations();
  
  void CleanUpExtrapolations();
  
  HelixClass * GetExtrapolationHelix(TrackExtended * track);
  
  // void PrintOutMerging(TrackExtended * firstTrackExt, TrackExtended * secondTrackExt, int iopt);
  
  // void GeneralSorting(int * index, float * val, int direct, int nVal);
  
  int SegmentRadialOverlap(TrackExtended* pTracki, TrackExtended* pTrackj);
  bool VetoMerge(TrackExtended* firstTrackExt, TrackExtended* secondTrackExt);
  
  
  int _nRun ;
  int _nEvt ;
  
  MarlinTrk::HelixFit* _fastfitter;
  
  /** pointer to the IMarlinTrkSystem instance 
   */
  MarlinTrk::IMarlinTrkSystem* _trksystem ;
  
  Gaudi::Property<bool> _MSOn{this, "MultipleScatteringOn", true};
  Gaudi::Property<bool> _ElossOn{this, "EnergyLossOn", true};
  Gaudi::Property<bool> _SmoothOn{this, "SmoothOn", true};
  
  TrackExtendedVec _allSiTracks;
  TrackExtendedVec _allTPCTracks;
  TrackExtendedVec _allCombinedTracks;
  TrackExtendedVec _allNonCombinedTPCTracks;
  TrackExtendedVec _allNonCombinedSiTracks;
  TrackExtendedVec _trkImplVec;
  TrackerHitExtendedVec _allTPCHits;
  TrackerHitExtendedVec _allVTXHits;
  TrackerHitExtendedVec _allFTDHits;
  TrackerHitExtendedVec _allSITHits;
  TrackerHitExtendedVec _allSETHits;
  TrackerHitExtendedVec _allETDHits;
  
  float PI, PIOVER2, TWOPI;
  
  float _bField;
  //float _chi2PrefitCut;
  Gaudi::Property<int> _debug{this, "Debug", 0};
  Gaudi::Property<int> _forceMerging{this, "ForceSiTPCMerging", 0};
  Gaudi::Property<int> _mergeTPCSegments{this, "ForceTPCSegmentsMerging", 0};
  Gaudi::Property<int> _cutOnTPCHits{this, "CutOnTPCHits", 35};
  Gaudi::Property<int> _cutOnSiHits{this, "CutOnSiHits", 4};
  Gaudi::Property<int> _assignVTXHits{this, "AssignVTXHits", 1};
  Gaudi::Property<int> _assignFTDHits{this, "AssignFTDHits", 1};
  Gaudi::Property<int> _assignSITHits{this, "AssignSITHits", 1};
  Gaudi::Property<int> _assignTPCHits{this, "AssignTPCHits", 1};
  Gaudi::Property<int> _assignSETHits{this, "AssignSETHits", 1};
  Gaudi::Property<int> _assignETDHits{this, "AssignETDHits", 1};
  Gaudi::Property<int> _forbidOverlapInZTPC{this, "ForbidOverlapInZTPC", 0};
  Gaudi::Property<int> _forbidOverlapInZComb{this, "ForbidOverlapInZComb", 0};

  //float _dPCutForMerging;
  Gaudi::Property<float> _d0CutForMerging{this, "D0CutForMerging", 500.};
  Gaudi::Property<float> _z0CutForMerging{this, "Z0CutForMerging", 1000.};
  Gaudi::Property<float> _dOmegaForMerging{this, "OmegaCutForMerging", 0.25};
  Gaudi::Property<float> _angleForMerging{this, "AngleCutForMerging", 0.10};
  Gaudi::Property<float> _chi2FitCut{this, "Chi2FitCut", 100.};
  Gaudi::Property<float> _d0CutForForcedMerging{this, "D0CutForForcedMerging", 50.};
  Gaudi::Property<float> _z0CutForForcedMerging{this, "Z0CutForForcedMerging", 200.};
  Gaudi::Property<float> _dOmegaForForcedMerging{this, "OmegaCutForForcedMerging", 0.15};
  Gaudi::Property<float> _angleForForcedMerging{this, "AngleCutForForcedMerging", 0.05};
  Gaudi::Property<float> _d0CutToMergeTPC{this, "D0CutToMergeTPCSegments", 100.};
  Gaudi::Property<float> _z0CutToMergeTPC{this, "Z0CutToMergeTPCSegments", 5000.};
  Gaudi::Property<float> _dPCutToMergeTPC{this, "DeltaPCutToMergeTPCSegments", 0.1};
  Gaudi::Property<float> _PtCutToMergeTPC{this, "PtCutToMergeTPCSegments", 1.2};
  Gaudi::Property<float> _cosThetaCutHighPtMerge{this, "cosThetaCutHighPtMerge", 0.99};
  Gaudi::Property<float> _cosThetaCutSoftHighPtMerge{this, "cosThetaCutSoftHighPtMerge", 0.998};
  Gaudi::Property<float> _momDiffCutHighPtMerge{this, "momDiffCutHighPtMerge", 0.01};
  Gaudi::Property<float> _momDiffCutSoftHighPtMerge{this, "momDiffCutSoftHighPtMerge", 0.25};
  Gaudi::Property<float> _hitDistanceCutHighPtMerge{this, "hitDistanceCutHighPtMerge", 25.};
  Gaudi::Property<float> _maxHitDistanceCutHighPtMerge{this, "maxHitDistanceCutHighPtMerge", 50.};
  Gaudi::Property<float> _maxFractionOfOutliersCutHighPtMerge{this, "maxFractionOfOutliersCutHighPtMerge", 0.95};
  Gaudi::Property<int>   _nHitsExtrapolation{this, "NHitsExtrapolation", 35};
  Gaudi::Property<float> _distCutForVTXHits{this, "VTXHitToTrackDistance", 1.5};
  Gaudi::Property<float> _distCutForFTDHits{this, "FTDHitToTrackDistance", 2.};
  Gaudi::Property<float> _distCutForSITHits{this, "SITHitToTrackDistance", 2.};
  Gaudi::Property<float> _distCutForSETHits{this, "SETHitToTrackDistance", 2.};
  Gaudi::Property<float> _distCutForETDHits{this, "ETDHitToTrackDistance", 10.};
  Gaudi::Property<float> _distCutForTPCHits{this, "TPCHitToTrackDistance", 15.};
  Gaudi::Property<float> _d0TrkCut{this, "CutOnTrackD0", 500.};
  Gaudi::Property<float> _z0TrkCut{this, "CutOnTrackZ0", 500.};
  Gaudi::Property<float> _initialTrackError_d0{this, "InitialTrackErrorD0",1e6};
  Gaudi::Property<float> _initialTrackError_phi0{this, "InitialTrackErrorPhi0",1e2};
  Gaudi::Property<float> _initialTrackError_omega{this, "InitialTrackErrorOmega",1e-4};
  Gaudi::Property<float> _initialTrackError_z0{this, "InitialTrackErrorZ0",1e6};
  Gaudi::Property<float> _initialTrackError_tanL{this, "InitialTrackErrorTanL",1e2};
  Gaudi::Property<float> _maxChi2PerHit{this, "MaxChi2PerHit", 1e2};
  Gaudi::Property<float> _minChi2ProbForSiliconTracks{this, "MinChi2ProbForSiliconTracks", 0.001};
  Gaudi::Property<float> _vetoMergeMomentumCut{this, "VetoMergeMomentumCut", 2.5};
  Gaudi::Property<float> _maxAllowedPercentageOfOutliersForTrackCombination{this, "MaxAllowedPercentageOfOutliersForTrackCombination", 0.3};
  Gaudi::Property<int>   _maxAllowedSiHitRejectionsForTrackCombination{this, "MaxAllowedSiHitRejectionsForTrackCombination", 2};

  //float _dPCutForForcedMerging;
  
  bool _runMarlinTrkDiagnostics = false;
  std::string _MarlinTrkDiagnosticsName;
  
  std::map<TrackExtended*,HelixClass*> _trackExtrapolatedHelix;
  std::set<TrackExtended*> _candidateCombinedTracks;
  
  UTIL::BitField64* _encoder;
  int getDetectorID(edm4hep::ConstTrackerHit hit) { _encoder->setValue(hit.getCellID()); return (*_encoder)[lcio::ILDCellID0::subdet]; }
  int getSideID(edm4hep::ConstTrackerHit hit)     { _encoder->setValue(hit.getCellID()); return (*_encoder)[lcio::ILDCellID0::side]; };
  int getLayerID(edm4hep::ConstTrackerHit hit)    { _encoder->setValue(hit.getCellID()); return (*_encoder)[lcio::ILDCellID0::layer]; };
  int getModuleID(edm4hep::ConstTrackerHit hit)   { _encoder->setValue(hit.getCellID()); return (*_encoder)[lcio::ILDCellID0::module]; };
  int getSensorID(edm4hep::ConstTrackerHit hit)   { _encoder->setValue(hit.getCellID()); return (*_encoder)[lcio::ILDCellID0::sensor]; };

  
  void setupGearGeom() ;
  
  double _tpc_inner_r;
  double _tpc_outer_r;
  double _tpc_pad_height;
  int    _tpc_nrows;
  
//  struct VXD_Layer {
//    int nLadders;
//    double phi0;
//    double dphi;
//    double senRMin;
//    double supRMin;
//    double length;
//    double width;
//    double offset;
//    double senThickness;
//    double supThickness;
//  };
//  std::vector<VXD_Layer> _VXDgeo;
  
  unsigned int _nLayersVTX;
  
//  struct SIT_Layer {
//    int nLadders;
//    double phi0;
//    double dphi;
//    double senRMin;
//    double supRMin;
//    double length;
//    double width;
//    double offset;
//    double senThickness;
//    double supThickness;
//  };
//  std::vector<SIT_Layer> _SITgeo;
  
  unsigned int _nLayersSIT;
  
  unsigned int _nLayersSET;

  
//  struct FTD_Disk {
//    int nPetals;
//    double phi0;
//    double dphi;
//    
//    double alpha;
//    double rInner;
//    double height;
//    double innerBaseLength;
//    double outerBaseLength;
//    double senThickness;
//    double supThickness;
//    
//    double senZPos_even_petal1;
//    double senZPos_even_petal2;
//    double senZPos_even_petal3;
//    double senZPos_even_petal4;
//    
//    double supZPos_even_petal1;
//    double supZPos_even_petal2;
//    double supZPos_even_petal3;
//    double supZPos_even_petal4;
//    
//    double senZPos_odd_petal1;
//    double senZPos_odd_petal2;
//    double senZPos_odd_petal3;
//    double senZPos_odd_petal4;
//    
//    double supZPos_odd_petal1;
//    double supZPos_odd_petal2;
//    double supZPos_odd_petal3;
//    double supZPos_odd_petal4;
//    
//    
//    
//  };
//  
//  std::vector<FTD_Disk> _FTDgeo;
  std::vector<float> _zLayerFTD;
  gear::GearMgr* gearMgr;

  unsigned int _nLayersFTD;
  int _nPhiFTD; 
  bool  _petalBasedFTDWithOverlaps;

  DataHandle<edm4hep::TrackerHitCollection> _TPCTrackerHitColHdl{"TPCTrackerHits", Gaudi::DataHandle::Reader, this};
  DataHandle<edm4hep::TrackerHitCollection> _FTDSpacePointColHdl{"FTDSpacePoints", Gaudi::DataHandle::Reader, this};
  DataHandle<edm4hep::TrackerHitCollection> _FTDPixelTrackerHitColHdl{"FTDPixelTrackerHits", Gaudi::DataHandle::Reader, this};
  DataHandle<edm4hep::TrackerHitCollection> _SITTrackerHitColHdl{"SITSpacePoints", Gaudi::DataHandle::Reader, this};
  DataHandle<edm4hep::TrackerHitCollection> _SETTrackerHitColHdl{"SETSpacePoints", Gaudi::DataHandle::Reader, this};
  DataHandle<edm4hep::TrackerHitCollection> _VTXTrackerHitColHdl{"VTXTrackerHits", Gaudi::DataHandle::Reader, this};
 
  DataHandle<edm4hep::TrackerHitCollection> _inSITRawColHdl{"SITTrackerHits", Gaudi::DataHandle::Reader, this};
  DataHandle<edm4hep::TrackerHitCollection> _inSETRawColHdl{"SETTrackerHits", Gaudi::DataHandle::Reader, this};
  DataHandle<edm4hep::TrackerHitCollection> _inFTDRawColHdl{"FTDStripTrackerHits", Gaudi::DataHandle::Reader, this};
  //DataHandle<edm4hep::SimTrackerHitCollection> _inVTXRawColHdl{"VXDCollection", Gaudi::DataHandle::Reader, this};

  DataHandle<edm4hep::TrackCollection> _TPCTrackColHdl{"ClupatraTracks", Gaudi::DataHandle::Reader, this};
  DataHandle<edm4hep::TrackCollection> _SiTrackColHdl{"SiTracks", Gaudi::DataHandle::Reader, this};

  DataHandle<edm4hep::TrackCollection> _OutputTrackColHdl{"MarlinTrkTracks", Gaudi::DataHandle::Writer, this};
};

#endif



