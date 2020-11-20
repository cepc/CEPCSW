// *********************************************************
// Implemented from Mokka
// *********************************************************
//
// $Id: TimeProjectionChamberSensitiveDetector.hh,v 1.1 2009/05/13 08:22:57 steve Exp $

#ifndef TimeProjectionChamberSensitiveDetector_h
#define TimeProjectionChamberSensitiveDetector_h

#include "DetSimSD/DDG4SensitiveDetector.h"
#include "DDG4/Defs.h"

class TimeProjectionChamberSensitiveDetector: public DDG4SensitiveDetector {
 public:
  TimeProjectionChamberSensitiveDetector(const std::string& name, dd4hep::Detector& description);
  
  void Initialize(G4HCofThisEvent* HCE);
  G4bool ProcessHits(G4Step* step, G4TouchableHistory* history);
  void EndOfEvent(G4HCofThisEvent* HCE);
  
  void setThreshold(double e){m_thresholdEnergyDeposit=e;};
  void setSameStepLimit(bool flag){m_lowPtStepLimit=(!flag);};
  void setWriteMCTruthForLowPtHits(bool flag){m_writeMCParticleForLowPtHits=flag;};
  void setLowPtCut(double cut){m_lowPtCut=cut;};
  void setLowPtMaxHitSeparation(double length){m_lowPtMaxHitSeparation=length;};
  
  /// helper function to avoid code duplication, writes a low Pt hit to the collection
  void DepositLowPtHit();
  
  /// helper function to avoid code duplication, resets all cumulative variables
  void ResetCumulativeVariables();
  
  /// helper function to avoid code duplication,
  /// adds energy, track length and momentum of a low pt step to the cumulative variables
  void CumulateLowPtStep(dd4hep::sim::Geant4StepHandler& h);
  
 protected:

  double m_thresholdEnergyDeposit = 0;
  bool   m_lowPtStepLimit = false;
  bool   m_writeMCParticleForLowPtHits = false;
  double m_lowPtCut = 0;
  double m_lowPtMaxHitSeparation = 0;
  bool   m_trackingPhysicsListELossOn = true;

  HitCollection* m_hc = nullptr;
  HitCollection* m_spaceHC = nullptr;
  HitCollection* m_lowPtHC = nullptr;
  
  dd4hep::Position      CrossingOfPadRingCentre;
  dd4hep::sim::Momentum MomentumAtPadRingCentre;
  double                dEInPadRow;
  double                globalTimeAtPadRingCentre;
  double                pathLengthInPadRow;
  double                CumulativePathLength;
  double                CumulativeEnergyDeposit;
  dd4hep::Position      CumulativeMeanPosition; 
  dd4hep::sim::Momentum CumulativeMeanMomentum; 
  int                   CumulativeNumSteps;
  
  dd4hep::Position PreviousPostStepPosition; //< the end point of the previous step
  int    CurrentPDGEncoding; //< the PDG encoding of the particle causing the cumulative energy deposit
  int    CurrentTrackID; //< the TrackID of the particle causing the cumulative energy deposit
  double CurrentGlobalTime; ///< the global time of the track causing the cumulative energy deposit
  int    CurrentCopyNumber; ///< copy number of the preStepPoint's TouchableHandle for the cumulative energy deposit
  
};
#endif
