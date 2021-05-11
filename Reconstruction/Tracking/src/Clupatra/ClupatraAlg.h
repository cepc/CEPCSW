#ifndef ClupatraProcessor_h
#define ClupatraProcessor_h 1

#include "assert.h"
#include "k4FWCore/DataHandle.h"
#include "GearSvc/IGearSvc.h"
#include "TrackSystemSvc/ITrackSystemSvc.h"

#include "GaudiAlg/GaudiAlgorithm.h"
#include "edm4hep/Track.h"
#include "edm4hep/TrackerHitConst.h"
#include "edm4hep/TrackCollection.h"
#include <string>

#include "gear/TPCModule.h"
#include "RuntimeMap.h"
#include "clupatra_new.h"
#include "Tracking/TrackingHelper.h"

#include "TrackSystemSvc/MarlinTrkUtils.h"
#include "TrackSystemSvc/HelixTrack.h"
#include "TrackSystemSvc/HelixFit.h"
#include "TrackSystemSvc/IMarlinTrack.h"

namespace MarlinTrk{
  class IMarlinTrkSystem ;
}

namespace gear{
  class TPCParameters ;
  class PadRowLayout2D ;
}

/** ClupatraAlgorithm : nearest neighbour clustering seeded pattern recognition using Kalman-Filter for extrapolation and
 *  hit search.
 *
 *   @parameter TPCHitCollection         Name of the tpc hit input collections
 *   @parameter OutputCollection         Name of the output collection with final TPC tracks
 *   @parameter SegmentCollectionName    Name of the output collection that has the individual track segments
 *   @parameter CreateDebugCollections   optionally create some debug collection with intermediate track segments and used and unused hits
 *
 *   @parameter DistanceCut              Cut for distance between hits in mm for the seed finding
 *   @parameter Chi2Cut                  the maximum chi2-distance for which a hit is considered for merging
 *   @parameter CosAlphaCut              Cut for max.angle between hits in consecutive layers for seed finding - NB value should be smaller than 1 - default is 0.9999999
 *   @parameter MaxDeltaChi2             the maximum delta Chi2 after filtering for which a hit is added to a track segement
 *
 *   @parameter DuplicatePadRowFraction  allowed fraction of hits in same pad row per track
 *   @parameter NLoopForSeeding          number of seed finding loops - every loop increases the distance cut by DistanceCut/NLoopForSeeding
 *   @parameter NumberOfZBins            number of bins in z over total length of TPC - hits from different z bins are nver merged
 *   @parameter PadRowRange              number of pad rows used in initial seed clustering
 *
 *   @parameter MaxStepWithoutHit                 the maximum number of layers without finding a hit before hit search search is stopped
 *   @parameter MinLayerFractionWithMultiplicity  minimum fraction of layers that have a given multiplicity, when forcing a cluster into sub clusters
 *   @parameter MinLayerNumberWithMultiplicity    minimum number of layers that have a given multiplicity, when forcing a cluster into sub clusters
 *   @parameter MinimumClusterSize                minimum number of hits per cluster
 *
 *   @parameter TrackEndsOuterCentralDist maximum radial distance [mm] from outer field cage of last hit, such that the track is considered to end at the end
 *   @parameter TrackEndsOuterForwardDist maximum distance in z [mm] from endplate of last hit, such that the track is considered to end at the end
 *   @parameter TrackIsCurlerOmega        minimum curvature omega of a track segment for being considered a curler
 *   @parameter TrackStartsInnerDist      maximum radial distance [mm] from inner field cage of first hit, such that the track is considered to start at the beginning
 *
 *   @parameter EnergyLossOn             Use Energy Loss in Fit
 *   @parameter MultipleScatteringOn     Use MultipleScattering in Fit
 *   @parameter SmoothOn                 Smooth All Mesurement Sites in Fit
 *
 *   @parameter pickUpSiHits             try to pick up hits from Si-trackers
 *   @parameter SITHitCollection         name of the SIT hit collections - used to extend TPC tracks if (pickUpSiHits==true)
 *   @parameter VXDHitCollection         name of the VXD hit collections - used to extend TPC tracks if (pickUpSiHits==true)
 *
 *   @parameter Verbosity               verbosity level of this processor ("DEBUG0-4,MESSAGE0-4,WARNING0-4,ERROR0-4,SILENT")
 *
 * @author F.Gaede, DESY, 2011/2012
 * @version $Id: ClupatraProcessor.h 4488 2013-04-05 12:03:59Z volynets $
 */
class ClupatraAlg : public GaudiAlgorithm {
    // friend class AlgFactory<ClupatraAlg>;
 public:


  ClupatraAlg(const std::string& name, ISvcLocator* svcLoc) ;

  /** Called at the begin of the job before anything is read.
   * Use to initialize the processor, e.g. book histograms.
   */
  virtual StatusCode initialize() ;

  /** Called for every run.
   */
  // virtual void processRunHeader( LCRunHeader* run ) ;

  /** Called for every event - the working horse.
   */
  virtual StatusCode execute() ;

  // virtual void check( LCEvent * evt ) ;

  /** Called after data processing for clean up.
   */
  virtual StatusCode finalize() ;


 protected:

  /** helper method to compute a few track segment parameters (start and end points, z spread,...)
   */
  void computeTrackInfo(edm4hep::ConstTrack lTrk) ;


  StatusCode pickUpSiTrackerHits(edm4hep::TrackCollection* trackCol) ;
  void printParameters();

  /** Input collection name.
   */
  std::string _colName ;
  std::string _vxdColName ;
  std::string _sitColName ;
  std::string _outColName ;
  std::string  _segmentsOutColName ;

  Gaudi::Property<float> _distCut{this, "DistanceCut", 40.0}; // Cut for distance between hits in mm for the seed finding
  Gaudi::Property<float> _cosAlphaCut{this, "CosAlphaCut", 0.9999999}; // Cut for max.angle between hits in consecutive layers for seed finding - NB value should be smaller than 1 - default is 0.9999999 !!!
  Gaudi::Property<int> _nLoop{this, "NLoopForSeeding", 4}; // number of seed finding loops - every loop increases the distance cut by DistanceCut/NLoopForSeeding
  Gaudi::Property<int> _minCluSize{this, "MinimumClusterSize", 6}; // minimum number of hits per cluster
  Gaudi::Property<float> _duplicatePadRowFraction{this, "DuplicatePadRowFraction", 0.1}; // allowed fraction of hits in same pad row per track
  Gaudi::Property<float> _dChi2Max{this, "MaxDeltaChi2", 35.}; // the maximum delta Chi2  after filtering for which a hit is added to a track segement
  Gaudi::Property<float> _chi2Cut{this, "Chi2Cut", 100.}; // the maximum chi2-distance for which a hit is considered for merging
  Gaudi::Property<int> _maxStep{this, "MaxStepWithoutHit", 6}; // the maximum number of layers without finding a hit before hit search search is stopped
  Gaudi::Property<bool> _pickUpSiHits{this, "pickUpSiHits", false};  // "try to pick up hits from Si-trackers"
  Gaudi::Property<bool> _createDebugCollections{this, "CreateDebugCollections", false}; // optionally create some debug collection with intermediate track segments and used and unused hits
  Gaudi::Property<int> _padRowRange{this, "PadRowRange", 15}; // "number of pad rows used in initial seed clustering"
  Gaudi::Property<int> _nZBins{this, "NumberOfZBins", 150}; // "number of bins in z over total length of TPC - hits from different z bins are nver merged"
  Gaudi::Property<float> _minLayerFractionWithMultiplicity{this, "MinLayerFractionWithMultiplicity" , 0.5}; // "minimum fraction of layers that have a given multiplicity, when forcing a cluster into sub clusters"
  Gaudi::Property<float> _minLayerNumberWithMultiplicity{this, "MinLayerNumberWithMultiplicity", 3}; // "minimum number of layers that have a given multiplicity, when forcing a cluster into sub clusters"
  Gaudi::Property<float> _trackStartsInnerDist{this, "TrackStartsInnerDist", (float) 25.}; // "maximum radial distance [mm] from inner field cage of first hit, such that the track is considered to start at the beginning "
  Gaudi::Property<float> _trackEndsOuterCentralDist{this, "TrackEndsOuterCentralDist", 25.}; // "maximum radial distance [mm] from outer field cage of last hit, such that the track is considered to end at the end " 
  Gaudi::Property<float> _trackEndsOuterForwardDist{this, "TrackEndsOuterForwardDist" , 40.}; // "maximum distance in z [mm] from endplate of last hit, such that the track is considered to end at the end " ,
  Gaudi::Property<float> _trackIsCurlerOmega{this, "TrackIsCurlerOmega", 0.001}; // "minimum curvature omega of a track segment for being considered a curler"
  Gaudi::Property<bool> _MSOn{this, "MultipleScatteringOn", false};
  Gaudi::Property<bool> _ElossOn{this, "EnergyLossOn", true};
  Gaudi::Property<bool> _SmoothOn{this, "SmoothOn", false};


  DataHandle<edm4hep::TrackerHitCollection> _TPCHitCollectionHandle{"TPCTrackerHits", Gaudi::DataHandle::Reader, this};
  DataHandle<edm4hep::TrackerHitCollection> _SITHitCollectionHandle{"SIDTrackerHits", Gaudi::DataHandle::Reader, this};
  DataHandle<edm4hep::TrackerHitCollection> _VTXHitCollectionHandle{"VTXTrackerHits", Gaudi::DataHandle::Reader, this};

  DataHandle<edm4hep::TrackCollection> _ClupatraTrackCollectionHandle{"ClupatraTracks", Gaudi::DataHandle::Writer, this};
  DataHandle<edm4hep::TrackCollection> _ClupatraTrackSegmentCollectionHandle{"ClupatraSegmentTracks", Gaudi::DataHandle::Writer, this};

  float _bfield ;

  int _nRun ;
  int _nEvt ;

  MarlinTrk::IMarlinTrkSystem* _trksystem ;

  const gear::TPCParameters*  _gearTPC ;

  TTree *tree;
  double omega;
  double pt;
  int totalCandidates;

//   NNClusterer* _clusterer ;

} ;

#endif



