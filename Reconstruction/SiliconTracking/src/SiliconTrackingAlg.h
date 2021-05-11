#ifndef SiliconTrackingAlg_h
#define SiliconTrackingAlg_h

//#include "marlin/Processor.h"
//#include <marlin/Global.h>
#include "k4FWCore/DataHandle.h"
#include "GaudiAlg/GaudiAlgorithm.h"
#include "edm4hep/EventHeaderCollection.h"
#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/SimTrackerHitCollection.h"
#include "edm4hep/TrackerHitCollection.h"
#include "edm4hep/TrackCollection.h"
#include "edm4hep/MCRecoTrackerAssociationCollection.h"

//#include "lcio.h"
#include <string>
#include <vector>
#include <cmath>
//#include <IMPL/TrackImpl.h>
#include "DataHelper/ClusterExtended.h"
#include "DataHelper/TrackExtended.h"
#include "DataHelper/TrackerHitExtended.h"
#include "DataHelper/HelixClass.h"

#include "TrackSystemSvc/IMarlinTrack.h"

#include <UTIL/BitField64.h>
#include <UTIL/ILDConf.h>

//using namespace edm4hep ;

namespace gear{
  class GearMgr ;
}

namespace MarlinTrk {
  class HelixFit;
  class IMarlinTrkSystem ;
}

namespace UTIL{
  class LCRelationNavigator ;
}

class TrackExtended;
class TrackerHitExtended;
class HelixClass;

/** === Silicon Tracking Processor === <br>
 * Processor performing stand-alone pattern recognition
 * in the vertex detector (VTX), forward tracking disks and SIT. <br>
 * The procedure consists of three steps : <br> 
 * 1) Tracking in VTX and SIT ; <br>
 * 2) Tracking in FTD ; <br>
 * 3) Merging compatible track segments reconstructed in VTX and FTD <br>
 * STEP 1 : TRACKING IN VTX and SIT <br>
 * Algorithm starts with finding of hit triplets satisfying helix hypothesis <br> 
 * in three different layers. Two layers of SIT are effectively considered as outermost <br>
 * layers of the vertex detector. To accelerate procedure, the 4-pi solid angle
 * is divided in NDivisionsInTheta and NDivisionsInPhi sectors in cosQ and Phi, 
 * respectively. Triplets are looked for in 2x2 window of adjacent sectors. 
 * Once triplet is found, attempt is made to associate additional hits to 
 * track. Combinatin of hits is accepted for further analysis if the Chi2 
 * of the fit is less than certain predefined threshold. All accepted 
 * combinations are sorted in ascending order of their Chi2. First track candidate 
 * in the sorted array is automatically accepted. The hits belonging to this track are 
 * marked as used, and track candidates sharing these hits are discarded.
 * The procedure proceeds with increasing index of track candidate in the sorted 
 * array until all track candidate have been output or discarded. <br>
 * STEP 2 : TRACKING IN FTD <br>
 * In the next step tracking in FTD is performed. The strategy of tracking in the FTD 
 * is the same as used for tracking in the VTX+SIT. <br>
 * STEP 3 : MERGING TRACK SEGMENTS FOUND IN FTD AND VTX+SIT <br>
 * In the last step, track segments reconstructed in the FTD and VTX+SIT, belonging to the
 * same track  are identified and merged into one track. All possible 
 * pairings are tested for their compatibility.
 * The number of pairings considered is Ntrk_VTX_SIT*Ntrk_FTD, where Ntrk_VTX_SIT is the number of 
 * track segments reconstructed in the first step in VTX+SIT (segments containing solely VTX and SIT hits) and
 * Ntrk_FTD is the number of track segments reconstructed in the second step 
 * (segments containing solely FTD hits).
 * Pair of segments is accepted for further examination if the angle between track segments and 
 * than certain specified threshold.
 * Pairing satisfying this condition is subjected for 
 * addtitional test. The fit is performed on unified array of hits belonging to both segments. 
 * If the chi2 of the fit does not exceed predefined cut value two segments are unified into 
 * one track. 
 * <h4>Input collections and prerequisites</h4> 
 * Processor requires collection of digitized vertex, sit and ftd tracker hits. <br>
 * If such a collections with the user specified names do not exist 
 * processor takes no action. <br>
 * <h4>Output</h4>
 * Processor produces an LCIO collection of the Tracks. Each track is characterised by 
 * five parameters : Omega (signed curvuture), Tan(lambda) where
 * lambda is the dip angle, Phi (azimuthal angle @ point of closest approach), D0 (signed impact parameter),
 * Z0 (displacement along z axis at the point of closest approach to IP). Covariance matrix for these parameters is also provided.
 * Only lower left corner of the covariance matrix is stored. The sequence of the covariance matrix elements 
 * assigned to track is the following: <br>
 * (Omega,Omega) <br>
 * (Omega,TanLambda), (TanLambda,TanLambda) <br>
 * (Omega,Phi), (TanLamda,Phi), (Phi,Phi) <br>
 * (Omega,D0), (TanLambda,D0), (Phi,D0), (D0,D0) <br>
 * (Omega,Z0), (TanLambda,Z0), (Phi,Z0), (D0,Z0), (Z0,Z0) <br>
 * The number of hits in the different subdetectors associated
 * with each track can be accessed via method Track::getSubdetectorHitNumbers().
 * This method returns vector of integers : <br>
 * number of VTX hits in track is the first element in this vector  
 * (Track::getSubdetectorHitNumbers()[0]) <br>
 * number of FTD hits in track is the second element in this vector  
 * (Track::getSubdetectorHitNumbers()[1]) <br>
 * number of SIT hits in track is the third element in this vector  
 * (Track::getSubdetectorHitNumbers()[2]) <br>
 * Output track collection has a name "SiTracks". <br>
 * @param VTXHitCollectionName name of input VTX TrackerHit collection <br>
 * (default parameter value : "VTXTrackerHits") <br>
 * @param FTDHitCollectionName name of input FTD TrackerHit collection <br>
 * (default parameter value : "FTDTrackerHits") <br>
 * @param SITHitCollectionName name of input SIT TrackerHit collection <br>
 * (default parameter value : "SITTrackerHits") <br>
 * @param SiTrackCollectionName name of the output Silicon track collection <br>
 * (default parameter value : "SiTracks") <br>
 * @param LayerCombinations combinations of layers used to search for hit triplets in VTX+SIT <br>
 * (default parameters : 6 4 3  6 4 2  6 3 2  5 4 3  5 4 2  5 3 2  4 3 2  4 3 1  4 2 1  3 2 1) <br> 
 * Note that in the VTX+SIT system the first and the second layers of SIT have indicies nLayerVTX and nLayerVTX+1. 
 * Combination given above means that triplets are looked first in layers 6 4 3, and then 
 * in 6 4 2;  5 4 3;  6 3 2 etc. NOTE THAT LAYER INDEXING STARTS FROM 0.
 * LAYER 0 is the innermost layer  <br>
 * @param LayerCombinationsFTD combinations of layers used to search for hit triplets in FTD <br>
 * (default parameters 6 5 4  5 4 3  5 4 2  5 4 1  5 3 2  5 3 1  5 2 1  4 3 2  4 3 1  
 *  4 3 0  4 2 1  4 2 0  4 1 0  3 2 1  3 2 0  3 1 0  2 1 0). 
 * NOTE THAT TRACKS IN FTD ARE SEARCHED ONLY IN ONE HEMISPHERE. TRACK IS NOT 
 * ALLOWED TO HAVE HITS BOTH IN BACKWARD AND FORWARD PARTS OF FTD SIMULTANEOUSLY. 
 * @param NDivisionsInPhi Number of divisions in Phi for tracking in VTX+SIT <br>
 * (default value is 40) <br>
 * @param NDivisionsInTheta Number of divisions in cosQ for tracking in VTX+SIT <br>
 * (default value is 40) <br>
 * @param NDivisionsInPhiFTD Number of divisions in Phi for tracking in FTD <br>
 * (default value is 3) <br>
 * @param Chi2WRphiTriplet weight on chi2 in R-Phi plane for track with 3 hits <br>
 * (default value is 1) <br>
 * @param Chi2WZTriplet weight on chi2 in S-Z plane for track with 3 hits <br>
 * (default value is 0.5) <br>
 * @param Chi2WRphiQuartet weight on chi2 in R-Phi plane to accept track with 4 hits <br>
 * (default value is 1) <br>
 * @param Chi2WZQuartet weight on chi2 in S-Z plane for track with 4 hits <br>
 * (default value is 0.5) <br>
 * @param Chi2WRphiSeptet weight on chi2 in R-Phi plane for track with 5 and more hits <br>
 * (default value is 1) <br>
 * @param Chi2WZSeptet Cut on chi2 in S-Z plane for track with 5 and more hits <br>
 * (default value is 0.5) <br>
 * @param Chi2FitCut Cut on chi2/ndf to accept track candidate <br>
 * (default value is 100.) <br>
 * @param AngleCutForMerging cut on the angle between two track segments.  
 * If the angle is greater than this cut, segments are not allowed to be merged. <br>
 * (default value is 0.1) <br>
 * @param MinDistCutAttach cut on the distance (in mm) from hit to the helix. This parameter is used
 * to decide whether hit can be attached to the track. If the distance is less than 
 * cut value. The track is refitted with a given hit being added to the list of hits already 
 * assigned for the track. Additional hit is assigned if chi2 of the new fit has good chi2. <br>
 * (default value is 2 ) <br>
 * @param MinLayerToAttach the minimal layer index to attach VTX hits to the found hit triplets <br>
 * (default value is -1) <br>
 * @param CutOnZ0 cut on Z0 parameter of track (in mm). If abs(Z0) is greater than the cut value, track is 
 * discarded (used to suppress fake
 * track rate in the presence of beam induced background hits) <br>
 * (default value is 100) <br>
 * @param CutOnD0 cut on D0 parameter of track (in mm). If abs(D0) is greater than the cut value, track is 
 * discarded (used to suppress fake
 * track rate in the presence of beam induced background hits) <br>
 * (default value is 100) <br>
 * @param CutOnPt cut on Pt (GeV/c). If Pt is less than this cut, track is discarded (used to suppress fake
 * track rate in the presence of beam induced background hits) <br>
 * (default value is 0.1) <br>
 * @param MinimalHits minimal number of hits in track required <br>
 * (default value is 3) <br>
 * @param NHitsChi2 Maximal number of hits for which a track with n hits is aways better than one with n-1 hits.
 * For tracks with equal or more than NHitsChi2 the track  with the lower \f$\chi^2\f$ is better.
 * (default value is 5) <br>
 * @param FastAttachment if this flag is set to 1, less accurate but fast procedure to merge additional hits to tracks is used <br> 
 * if set to 0, a more accurate, but slower procedure is invoked <br>
 * (default value is 0) <br>
 * @param UseSIT When this flag is set to 1, SIT is included in pattern recognition. When this flag is set
 * to 0, SIT is excluded from the procedure of pattern recognition <br>
 * (default value is 1) <br>
 * <br>
 * @author A. Raspereza (MPI Munich)<br>
 */
class SiliconTrackingAlg : public GaudiAlgorithm {
 public:
  
  SiliconTrackingAlg(const std::string& name, ISvcLocator* svcLoc);
  
  virtual StatusCode initialize() ;
  
  virtual StatusCode execute() ; 
  
  virtual StatusCode finalize() ;
  
 protected:
  
  int _nRun ;
  int _nEvt ;
  //EVENT::LCEvent* _current_event;
  int _nLayers;
  unsigned int _nLayersVTX;
  unsigned int _nLayersSIT;
  int _ntriplets, _ntriplets_good, _ntriplets_2MCP, _ntriplets_3MCP, _ntriplets_1MCP_Bad, _ntriplets_bad;
  
  MarlinTrk::HelixFit* _fastfitter;
  gear::GearMgr* _GEAR;
  /** pointer to the IMarlinTrkSystem instance 
   */
  MarlinTrk::IMarlinTrkSystem* _trksystem ;
  //bool _runMarlinTrkDiagnostics;
  //std::string _MarlinTrkDiagnosticsName;
  typedef std::vector<int> IntVec;
  
  Gaudi::Property<IntVec> _Combinations{this, "LayerCombinations",
      {8,6,5, 8,6,4, 8,6,3, 8,6,2, 8,5,3, 8,5,2, 8,4,3, 8,4,2, 6,5,3, 6,5,2, 6,4,3, 6,4,2, 6,3,1, 6,3,0, 6,2,1, 6,2,0,
	5,3,1, 5,3,0, 5,2,1, 5,2,0, 4,3,1, 4,3,0, 4,2,1, 4,2,0}};
  Gaudi::Property<IntVec> _CombinationsFTD{this, "LayerCombinationsFTD", {4,3,2, 4,3,1, 4,3,0, 4,2,1, 4,2,0, 4,1,0, 3,2,1, 3,2,0, 3,1,0, 2,1,0,
	9,8,7, 9,8,6, 9,8,5, 9,7,6, 9,7,5, 9,6,5, 8,7,6, 8,7,5, 8,6,5, 7,6,5}};
  Gaudi::Property<int> _nDivisionsInPhi{this, "NDivisionsInPhi", 80};
  Gaudi::Property<int> _nDivisionsInPhiFTD{this, "NDivisionsInPhiFTD", 30};
  Gaudi::Property<int> _nDivisionsInTheta{this, "NDivisionsInTheta", 80};
  Gaudi::Property<float> _chi2WRPhiTriplet{this, "Chi2WRphiTriplet", 1.};
  Gaudi::Property<float> _chi2WRPhiQuartet{this, "Chi2WRphiQuartet", 1.};
  Gaudi::Property<float> _chi2WRPhiSeptet{this, "Chi2WRphiSeptet", 1.};
  Gaudi::Property<float> _chi2WZTriplet{this, "Chi2WZTriplet", 0.5};
  Gaudi::Property<float> _chi2WZQuartet{this, "Chi2WZQuartet", 0.5};
  Gaudi::Property<float> _chi2WZSeptet{this, "Chi2WZSeptet", 0.5};
  Gaudi::Property<float> _chi2FitCut{this, "Chi2FitCut", 120.};
  Gaudi::Property<float> _angleCutForMerging{this, "AngleCutForMerging", 0.1};
  Gaudi::Property<float> _minDistCutAttach{this, "MinDistCutAttach", 2.5};
  Gaudi::Property<float> _minimalLayerToAttach{this, "MinLayerToAttach", -1};
  Gaudi::Property<float> _cutOnD0{this, "CutOnD0", 100.0};
  Gaudi::Property<float> _cutOnZ0{this, "CutOnZ0", 100.0};
  Gaudi::Property<float> _cutOnPt{this, "CutOnPt", 0.05};
  Gaudi::Property<int> _minimalHits{this, "MinimalHits",3};
  Gaudi::Property<int> _nHitsChi2{this, "NHitsChi2", 5};
  Gaudi::Property<int> _max_hits_per_sector{this, "MaxHitsPerSector", 100};
  Gaudi::Property<int> _attachFast{this, "FastAttachment", 0};
  Gaudi::Property<bool> _useSIT{this, "UseSIT", true};
  Gaudi::Property<float> _initialTrackError_d0{this, "InitialTrackErrorD0",1e6};
  Gaudi::Property<float> _initialTrackError_phi0{this, "InitialTrackErrorPhi0",1e2};
  Gaudi::Property<float> _initialTrackError_omega{this, "InitialTrackErrorOmega",1e-4};
  Gaudi::Property<float> _initialTrackError_z0{this, "InitialTrackErrorZ0",1e6};
  Gaudi::Property<float> _initialTrackError_tanL{this, "InitialTrackErrorTanL",1e2};
  Gaudi::Property<float> _maxChi2PerHit{this, "MaxChi2PerHit", 1e2};
  Gaudi::Property<int> _checkForDelta{this, "CheckForDelta", 1};
  Gaudi::Property<float> _minDistToDelta{this, "MinDistToDelta", 0.25};
  Gaudi::Property<bool> _MSOn{this, "MultipleScatteringOn", true};
  Gaudi::Property<bool> _ElossOn{this, "EnergyLossOn", true};
  Gaudi::Property<bool> _SmoothOn{this, "SmoothOn", true};
  Gaudi::Property<float> _helix_max_r{this, "HelixMaxR", 2000.};
  
  //std::vector<int> _colours;  
  
  /** helper function to get collection using try catch block */
  //LCCollection* GetCollection(  LCEvent * evt, std::string colName ) ;
  
  /** helper function to get relations using try catch block */
  //LCRelationNavigator* GetRelations( LCEvent * evt, std::string RelName ) ;
  
  /** input MCParticle collection and threshold used for Drawing
   */
  //Gaudi::Property<Float> _MCpThreshold{this, "MCpThreshold", 0.1};
  //std::string  _colNameMCParticles;
  
  /// Compare tracks according to their chi2/ndf
  struct compare_TrackExtended{
    // n.b.: a and b should be TrackExtended const *, but the getters are not const :-(
    bool operator()(TrackExtended *a, TrackExtended *b) const {
      if ( a == b ) return false;
      return (a->getChi2()/a->getNDF() < b->getChi2()/b->getNDF() );
    }
  };
  
  
  //std::string _VTXHitCollection;
  //std::string _FTDPixelHitCollection;
  //std::string _FTDSpacePointCollection;
  //std::string _SITHitCollection;
  //std::string _siTrkCollection;
  
  //std::vector< LCCollection* > _colTrackerHits;
  //std::map< LCCollection*, std::string > _colNamesTrackerHits;

  // Input collections
  DataHandle<edm4hep::EventHeaderCollection> _headerColHdl{"EventHeaderCol", Gaudi::DataHandle::Reader, this};
  DataHandle<edm4hep::MCParticleCollection> _inMCColHdl{"MCParticle", Gaudi::DataHandle::Reader, this};
  DataHandle<edm4hep::TrackerHitCollection> _inVTXColHdl{"VXDTrackerHits", Gaudi::DataHandle::Reader, this};
  DataHandle<edm4hep::TrackerHitCollection> _inFTDPixelColHdl{"FTDPixelTrackerHits", Gaudi::DataHandle::Reader, this};
  DataHandle<edm4hep::TrackerHitCollection> _inFTDSpacePointColHdl{"FTDSpacePoints", Gaudi::DataHandle::Reader, this};
  DataHandle<edm4hep::TrackerHitCollection> _inSITColHdl{"SITSpacePoints", Gaudi::DataHandle::Reader, this};
  DataHandle<edm4hep::TrackerHitCollection> _inSITRawColHdl{"SITTrackerHits", Gaudi::DataHandle::Reader, this};
  DataHandle<edm4hep::TrackerHitCollection> _inFTDRawColHdl{"FTDStripTrackerHits", Gaudi::DataHandle::Reader, this};
  // Output collections
  DataHandle<edm4hep::TrackCollection> _outColHdl{"SiTracks", Gaudi::DataHandle::Writer, this};
  //DataHandle<edm4hep::LCRelationCollection> _outRelColHdl{"TrackerHitRelations", Gaudi::DataHandle::Reader, this};
  
  std::vector<TrackerHitExtendedVec> _sectors;
  std::vector<TrackerHitExtendedVec> _sectorsFTD;
  
  /**
   * A helper class to allow good code readability by accessing tracks with N hits.
   * As the smalest valid track contains three hits, but the first index in a vector is 0,
   * this class hides the index-3 calculation. As the vector access is inline there should be
   * no performance penalty.
   */
  class TracksWithNHitsContainer {
  public:
    /// Empty all the vectors and delete the tracks contained in it.
    void clear();
    
    /// Set the size to allow a maximum of maxHit hits.
    inline void resize(size_t maxHits) {
      _tracksNHits.resize(maxHits-2);
      _maxIndex=(maxHits-3);
    }
    
    // Sort all track vectors according to chi2/nDof
    //      void sort();
    
    /// Returns the  TrackExtendedVec for track with n hits. 
    /// In case n is larger than the maximal number the vector with the largest n ist returned.
    /// \attention The smallest valid number is three! For
    /// performance reasons there is no safety check!
    inline TrackExtendedVec & getTracksWithNHitsVec( size_t nHits ) {
      //return _tracksNHits[ std::min(nHits-3, _maxIndex) ];
      // for debugging: with boundary check
      return _tracksNHits.at(std::min(nHits-3, _maxIndex));
    }
    
  protected:
    std::vector< TrackExtendedVec > _tracksNHits;
    size_t _maxIndex; /// local cache variable to avoid calculation overhead
  };
  
  TracksWithNHitsContainer _tracksWithNHitsContainer;
  
  int InitialiseVTX();
  int InitialiseFTD();
  void ProcessOneSector(int iSectorPhi, int iSectorTheta);
  void CleanUp();
  TrackExtended * TestTriplet(TrackerHitExtended * outerHit, 
                              TrackerHitExtended * middleHit,
                              TrackerHitExtended * innerHit,
                              HelixClass & helix);
  
  int BuildTrack(TrackerHitExtended * outerHit, 
                 TrackerHitExtended * middleHit,
                 TrackerHitExtended * innerHit,
                 HelixClass & helix, 
                 int innerlayer,
                 int iPhiLow, int iPhiUp,
                 int iTheta, int iThetaUp,
                 TrackExtended * trackAR);
  
  void Sorting( TrackExtendedVec & trackVec);
  void CreateTrack(TrackExtended * trackAR );
  void AttachRemainingVTXHitsSlow();
  void AttachRemainingFTDHitsSlow();
  void AttachRemainingVTXHitsFast();
  void AttachRemainingFTDHitsFast();
  void TrackingInFTD();
  int BuildTrackFTD(TrackExtended* trackAR, int* nLR, int iS);
  int AttachHitToTrack(TrackExtended * trackAR, TrackerHitExtended * hit, int iopt);
  
  void FinalRefit(edm4hep::TrackCollection*);//, edm4hep::LCRelationCollection*);
  
  float _bField;
  
  // two pi is not a constant in cmath. Calculate it, once!
  static const double TWOPI;
  
  double _dPhi;
  double _dTheta;
  double _dPhiFTD;
  
  float _resolutionRPhiVTX;
  float _resolutionZVTX;
  
  float _resolutionRPhiFTD;
  float _resolutionZFTD;
  
  float _resolutionRPhiSIT;
  float _resolutionZSIT;
  
  float _phiCutForMerging;
  float _tanlambdaCutForMerging;
  //float _angleCutForMerging;
  
  float _distRPhi;
  float _distZ;
  
  float _cutOnOmega;

  TrackExtendedVec _trackImplVec;
      
  int _nTotalVTXHits,_nTotalFTDHits,_nTotalSITHits;
  
  //  int _createMap;
  
  UTIL::BitField64* _encoder;
  int getDetectorID(edm4hep::ConstTrackerHit hit) { _encoder->setValue(hit.getCellID()); return (*_encoder)[lcio::ILDCellID0::subdet]; }
  int getSideID(edm4hep::ConstTrackerHit hit)     { _encoder->setValue(hit.getCellID()); return (*_encoder)[lcio::ILDCellID0::side]; };
  int getLayerID(edm4hep::ConstTrackerHit hit)    { _encoder->setValue(hit.getCellID()); return (*_encoder)[lcio::ILDCellID0::layer]; };
  int getModuleID(edm4hep::ConstTrackerHit hit)   { _encoder->setValue(hit.getCellID()); return (*_encoder)[lcio::ILDCellID0::module]; };
  int getSensorID(edm4hep::ConstTrackerHit hit)   { _encoder->setValue(hit.getCellID()); return (*_encoder)[lcio::ILDCellID0::sensor]; };
  
  StatusCode setupGearGeom() ;
  
  std::vector<float> _zLayerFTD;
  
  unsigned int _nlayersFTD;
  bool _petalBasedFTDWithOverlaps;

  int _output_track_col_quality;
  static const int _output_track_col_quality_GOOD;
  static const int _output_track_col_quality_FAIR;
  static const int _output_track_col_quality_POOR;
} ;

#endif



