#ifndef TrackSubsetAlg_h
#define TrackSubsetAlg_h 1

#include "k4FWCore/DataHandle.h"
#include "GaudiAlg/GaudiAlgorithm.h"

#include "edm4hep/TrackCollection.h"
//#include "edm4hep/Track.h"
#include "TrackSystemSvc/IMarlinTrkSystem.h"

#include "Math/ProbFunc.h"

/**  Processor that takes tracks from multiple sources and outputs them (or modified versions, or a subset of them)
 * as one track collection.
 * 
 *  <h4>Input - Prerequisites</h4>
 *  Track collections
 *
 *  <h4>Output</h4> 
 *  A single track collection
 * 
 * @param TrackInputCollections A vector of the input track collections <br>
 * (default value ForwardTracks SiTracks )
 * 
 * @param TrackOutputCollection Name of the output track collection <br>
 * (default value SubsetTracks )
 * 
 * @param MultipleScatteringOn Whether to take multiple scattering into account when fitting the tracks<br>
 * (default value true )
 *  
 * @param EnergyLossOn Whether to take energyloss into account when fitting the tracks<br>
 * (default value true )
 * 
 * @param SmoothOn Whether to smooth all measurement sites in fit<br>
 * (default value false )
 * 
 * @param Omega The parameter omega for the HNN. Controls the influence of the quality indicator. Between 0 and 1:
 * 1 means high influence of quality indicator, 0 means no influence. 
 * 
 * @author Robin Glattauer, HEPHY
 * 
 */

class TrackSubsetAlg : public GaudiAlgorithm {
 public:

  TrackSubsetAlg(const std::string& name, ISvcLocator* svcLoc);
  
  virtual StatusCode initialize() ;

  virtual StatusCode execute() ;

  virtual StatusCode finalize() ;
  
 protected:
  MarlinTrk::IMarlinTrkSystem* _trkSystem;
  /* Input collection */
  std::vector<DataHandle<edm4hep::TrackCollection>* > _inTrackColHdls;
  std::vector<DataHandle<edm4hep::TrackerHitCollection>* > _inTrackerHitColHdls;
  /* Output collection */
  DataHandle<edm4hep::TrackCollection> _outColHdl{"SubsetTracks", Gaudi::DataHandle::Writer, this};
  
  Gaudi::Property<std::vector<std::string> > _trackInputColNames{this, "TrackInputCollections", {"ForwardTracks", "SiTracks"}};
  Gaudi::Property<std::vector<std::string> > _trackerHitInputColNames{this, "RawTrackerHitCollections",
      {"VXDTrackerHits", "SITTrackerHits", "FTDPixelTrackerHits", "FTDStripTrackerHits"}};
  Gaudi::Property<bool> _MSOn{this, "MultipleScatteringOn", true};
  Gaudi::Property<bool> _ElossOn{this, "EnergyLossOn", true};
  Gaudi::Property<bool> _SmoothOn{this, "SmoothOn", false};
  Gaudi::Property<float> _initialTrackError_d0{this, "InitialTrackErrorD0",1e6};
  Gaudi::Property<float> _initialTrackError_phi0{this, "InitialTrackErrorPhi0",1e2};
  Gaudi::Property<float> _initialTrackError_omega{this, "InitialTrackErrorOmega",1e-4};
  Gaudi::Property<float> _initialTrackError_z0{this, "InitialTrackErrorZ0",1e6};
  Gaudi::Property<float> _initialTrackError_tanL{this, "InitialTrackErrorTanL",1e2};
  Gaudi::Property<double> _maxChi2PerHit{this, "MaxChi2PerHit", 1e2};
  Gaudi::Property<double> _omega{this, "Omega", 0.75};
  
  float _bField;
  
  int _nRun ;
  int _nEvt ;
};

/** A functor to return whether two tracks are compatible: The criterion is if the share a TrackerHit or more */
class TrackCompatibility{
 public:
  inline bool operator()( edm4hep::Track* trackA, edm4hep::Track* trackB ){
    unsigned nHitsA = trackA->trackerHits_size();
    unsigned nHitsB = trackB->trackerHits_size();
    for( unsigned i=0; i < nHitsA; i++){
      for( unsigned j=0; j < nHitsB; j++){
	if ( trackA->getTrackerHits(i) == trackB->getTrackerHits(j) ) return false;      // a hit is shared -> incompatible
      }
    }
    
    return true;      
  }
};


/** A functor to return the quality of a track, which is currently the chi2 probability. */
class TrackQI{
 public:
  /** @param trkSystem a pointer to an IMarlinTrkSystem, needed for fitting of tracks */
  TrackQI( MarlinTrk::IMarlinTrkSystem* trkSystem ): _trkSystem(trkSystem){}
  
  inline double operator()( edm4hep::Track* track ){
    return ROOT::Math::chisquared_cdf_c( track->getChi2() , track->getNdf() );   
  }
  
 protected:
  
  MarlinTrk::IMarlinTrkSystem* _trkSystem;
};
#endif



