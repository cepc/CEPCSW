#include "TimeProjectionChamberSensitiveDetector.h"

#include "G4Step.hh"
#include "G4VProcess.hh"
//#include "G4SteppingManager.hh"
#include "G4SDManager.hh"
//#include "UserTrackInformation.hh"
#include "DD4hep/DD4hepUnits.h"

TimeProjectionChamberSensitiveDetector::TimeProjectionChamberSensitiveDetector(const std::string& name,
									       dd4hep::Detector& description)
  : DDG4SensitiveDetector(name, description){
  
  G4String CollName1=name+"Collection";
  collectionName.insert(CollName1);
  
  G4String CollName2=name+"SpacePointCollection";
  collectionName.insert(CollName2);
  
  G4String CollName3=name+"LowPtCollection";
  collectionName.insert(CollName3);
}

void TimeProjectionChamberSensitiveDetector::Initialize(G4HCofThisEvent* HCE){
  m_hc = new HitCollection(GetName(), collectionName[0]);
  int HCID = G4SDManager::GetSDMpointer()->GetCollectionID(m_hc);
  HCE->AddHitsCollection( HCID, m_hc ); 

  m_spaceHC = new HitCollection(GetName(), collectionName[1]);
  int HCIDSpace = G4SDManager::GetSDMpointer()->GetCollectionID(m_spaceHC);
  HCE->AddHitsCollection( HCIDSpace, m_spaceHC );

  m_lowPtHC = new HitCollection(GetName(), collectionName[2]);
  int HCIDLowPt = G4SDManager::GetSDMpointer()->GetCollectionID(m_lowPtHC);
  HCE->AddHitsCollection( HCIDLowPt, m_lowPtHC );
  
  CumulativeNumSteps=0;
  
  dEInPadRow = 0.0;
  CumulativeEnergyDeposit = 0.0;
  globalTimeAtPadRingCentre=0.0;
  pathLengthInPadRow=0.0;
  CumulativePathLength=0.0;
  CrossingOfPadRingCentre.SetXYZ(0,0,0);// = dd4hep::Position(0,0,0);
  MomentumAtPadRingCentre.SetXYZ(0,0,0);// = dd4hep::sim::Momentum(0,0,0);
  
  CumulativeMeanPosition.SetXYZ(0,0,0);// = dd4hep::Position(0,0,0);
  CumulativeMeanMomentum.SetXYZ(0,0,0);// = dd4hep::sim::Momentum(0,0,0);
  
  PreviousPostStepPosition.SetXYZ(0,0,0);// = dd4hep::Position(0,0,0);
  CurrentPDGEncoding=0;
  CurrentTrackID=-1;
  CurrentGlobalTime=0.0;
  CurrentCopyNumber=0;
  
  /*
  if( m_lowPtStepLimit ) {
    std::cout << "TimeProjectionChamberSensitiveDetector: TPCLowPtStepLimit ON" << std::endl;
  }
  else {
    std::cout << "TimeProjectionChamberSensitiveDetector: TPCLowPtStepLimit OFF" << std::endl;
  }

  if( m_writeMCParticleForLowPtHits ) {
    std::cout << "TimeProjectionChamberSensitiveDetector: TPCLowPtStoreMCPForHits ON" << std::endl;
  }
  else {
    std::cout << "TimeProjectionChamberSensitiveDetector: TPCLowPtStoreMCPForHits OFF" << std::endl;
  }

  std::cout << "TimeProjectionChamberSensitiveDetector: Threshold for Energy Deposit = " << m_thresholdEnergyDeposit / CLHEP::eV << " eV" << G4endl;
  std::cout << "TimeProjectionChamberSensitiveDetector: TPCLowPtCut = " << m_lowPtCut / CLHEP::MeV << " MeV "<< G4endl;
  std::cout << "TimeProjectionChamberSensitiveDetector: TPCLowPt Max hit separation "<< m_lowPtMaxHitSeparation / CLHEP::mm << " mm" << G4endl;
  */
}

G4bool TimeProjectionChamberSensitiveDetector::ProcessHits(G4Step* step, G4TouchableHistory*){
  // FIXME: 
  // (i) in the following algorithm if a particle "appears" within a pad-ring half and 
  // leaves without passing through the middle of the pad-ring it will not create a hit in 
  // this ring
  // (ii) a particle that crosses the boundry between two pad-ring halves will have the hit 
  // placed on this surface at the last crossing point, and will be assinged the total energy 
  // deposited in the whole pad-ring. This is a possible source of bias for the hit
  
  G4TouchableHandle touchPost = step->GetPostStepPoint()->GetTouchableHandle(); 
  G4TouchableHandle touchPre  = step->GetPreStepPoint()->GetTouchableHandle(); 
  dd4hep::sim::Geant4StepHandler h(step);
  
  if (fabs(h.trackDef()->GetPDGCharge()) < 0.01) return true;
  
  const dd4hep::Position PrePosition   = h.prePos();
  const dd4hep::Position PostPosition  = h.postPos();
  const dd4hep::sim::Momentum Momentum = h.postMom(); //use post step momentum, Mokka history
  
  float ptSQRD = Momentum.X()*Momentum.X()+Momentum.Y()*Momentum.Y();
  if( ptSQRD >= (m_lowPtCut*m_lowPtCut) ){
    // Step finishes at a geometric boundry
    if(step->GetPostStepPoint()->GetStepStatus() == fGeomBoundary){
    //if(h.postStepStatus() == "GeomBoundary"){
      // step within the same pair of upper and lower pad ring halves
      if(int(touchPre->GetCopyNumber()/2)==int(touchPost->GetCopyNumber()/2)) {
        //this step must have ended on the boundry between these two pad ring halfs 
        //record the tracks coordinates at this position 
        //and return

        CrossingOfPadRingCentre.SetXYZ(PostPosition.X(),PostPosition.Y(),PostPosition.Z());
        MomentumAtPadRingCentre   = Momentum;  
        dEInPadRow               += h.deposit();
        globalTimeAtPadRingCentre = h.trkTime();
        pathLengthInPadRow       += h.stepLength();
        
        //	    G4cout << "step must have ended on the boundry between these two pad ring halfs" << G4endl;
        //	    G4cout << "CrossingOfPadRingCentre = "   
        //		   << sqrt( CrossingOfPadRingCentre[0]*CrossingOfPadRingCentre[0]
        //			    +CrossingOfPadRingCentre[1]*CrossingOfPadRingCentre[1] )
        //		   << G4endl;
        
        return true;
      }
      else if(!(CrossingOfPadRingCentre.X()==0.0 && CrossingOfPadRingCentre.Y()==0.0 && CrossingOfPadRingCentre.Z()==0.0)) {
        // the above IF statment is to catch the case where the particle "appears" in this pad-row half volume and 
        // leaves with out crossing the pad-ring centre, as mentioned above
        
        //it is leaving the pad ring couplet
        //write out a hit
        //make sure particle is added to MC list
        //and return
        //	    G4cout << "step must be leaving the pad ring couplet" << G4endl;
        //	    G4cout << "write out hit at " 
        //		   << sqrt( CrossingOfPadRingCentre[0]*CrossingOfPadRingCentre[0]
        //			    +CrossingOfPadRingCentre[1]*CrossingOfPadRingCentre[1] )
        //		   << " " << "dEdx = " << step->GetTotalEnergyDeposit()+dEInPadRow 
        //		   << " " << "step length = " << step->GetStepLength()+pathLengthInPadRow  
        //		   << G4endl;
        
        double dE = h.deposit()+dEInPadRow;
	if ( dE > m_thresholdEnergyDeposit || m_trackingPhysicsListELossOn == false ) {
	  // needed for delta electrons: all particles causing hits have to be saved
	  // TODO: set track flag to save

	  dd4hep::sim::Geant4TrackerHit* hit = new dd4hep::sim::Geant4TrackerHit(h.trkID(),h.trkPdgID(),dE,globalTimeAtPadRingCentre);
	  hit->cellID        = touchPre->GetCopyNumber()/2+1;//getCellID(step);
	  hit->energyDeposit = dE;
	  hit->position      = CrossingOfPadRingCentre;
	  hit->momentum      = MomentumAtPadRingCentre;
	  hit->length        = h.stepLength()+pathLengthInPadRow;
          m_hc->insert(hit);
        }
        
        // zero cumulative variables 
        dEInPadRow                = 0.0;
        globalTimeAtPadRingCentre = 0.0;
        pathLengthInPadRow        = 0.0;
        CrossingOfPadRingCentre.SetXYZ(0.,0.,0.);
        MomentumAtPadRingCentre.SetXYZ(0.,0.,0.);
        return true;
      }
    }
    //case for which the step remains within geometric volume
    //FIXME: need and another IF case to catch particles which Stop within the padring
    else if(step->GetPostStepPoint()->GetStepStatus() != fGeomBoundary){
    //else if(h.postStepStatus() != "GeomBoundary") {
      //the step is not limited by the step length 
      if(step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()!="StepLimiter"){
        // if(particle not stoped){ 	    
        //add the dEdx and return 
        //	    G4cout << "Step ended by Physics Process: Add dEdx and carry on" << G4endl;
        dEInPadRow         += h.deposit();
        pathLengthInPadRow += h.stepLength();
	return true;
        //}
        //else{
        //  write out the hit and clear counters
        //}
      }
      else {
	//double position_x = PostPosition[0];
        //double position_y = PostPosition[1];
        //double position_z = PostPosition[2];
        
        const dd4hep::sim::Momentum PreMomentum = h.preMom();
        const dd4hep::sim::Momentum PostMomentum = h.postMom();
        
        //double momentum_x = PostMomentum[0];
        //double momentum_y = PostMomentum[1];
        //double momentum_z = PostMomentum[2];
        
        
        //	    G4cout << "step must have been stopped by the step limiter" << G4endl;
        //	    G4cout << "write out hit at " 
        //		   << sqrt( position_x*position_x
        //			    +position_y*position_y )
        //		   << " " << "dEdx = " << step->GetTotalEnergyDeposit() 
        //		   << " " << "step length = " << step->GetStepLength()  
        //		   << G4endl;
        
        // write out step limited hit 
        // these are just space point hits so do not save the dE, which is set to ZERO
        //	    if ( step->GetTotalEnergyDeposit() > fThresholdEnergyDeposit ) 

        // needed for delta electrons: all particles causing hits have to be saved in the LCIO file
	// TODO: set track flag to be saved
        //UserTrackInformation *info = (UserTrackInformation*)(step->GetTrack()->GetUserInformation());
        //if (info) info->GetTheTrackSummary()->SetToBeSaved();

	dd4hep::sim::Geant4TrackerHit* hit = new dd4hep::sim::Geant4TrackerHit(h.trkID(),h.trkPdgID(),0.0,h.trkTime());
	hit->cellID        = touchPre->GetCopyNumber()/2+1;//getCellID(step);
	hit->energyDeposit = 0.0;
	hit->position      = PostPosition;
	hit->momentum      = PostMomentum;
	hit->length        = h.stepLength();
	m_spaceHC->insert(hit);
	
        // add dE and pathlegth and return
        dEInPadRow += h.deposit();
        pathLengthInPadRow += h.stepLength();
        return true;
      }
    }
  }
  else if(m_lowPtStepLimit) { // low pt tracks will be treated differently as their step length is limited by the special low pt steplimiter
    if ( ( PreviousPostStepPosition - PrePosition ).R() > 1.0e-6 * dd4hep::mm ) {
      // This step does not continue the previous path. Deposit the energy and begin a new Pt hit.
      if (CumulativeEnergyDeposit > m_thresholdEnergyDeposit) {
        DepositLowPtHit();
      }
      else {
        // reset the cumulative variables if the hit has not been deposited.
        // The previous track has ended and the cumulated energy left at the end 
        // was not enough to ionize
        //G4cout << "reset due to new track , discarding " << CumulativeEnergyDeposit / eV << " eV" << std::endl;
        ResetCumulativeVariables();
      }
    }
    else {
      //G4cout << "continuing track" << endl;
    }
    
    CumulateLowPtStep(h);  

    // check whether to deposit the hit
    if( ( CumulativePathLength > m_lowPtMaxHitSeparation )  ) {
      // hit is deposited because the step limit is reached and there is enough energy
      // to ionize
      if ( CumulativeEnergyDeposit > m_thresholdEnergyDeposit) {
        DepositLowPtHit();
      }
      //else {
      //G4cout << "not deposited, energy is " << CumulativeEnergyDeposit/eV << " eV" << std::endl;
      //}
    }
    else { // even if the step lenth has not been reached the hit might
           // be deposited because the particle track ends
      if ( h.post->GetKineticEnergy() == 0 ) {
	// only deposit the hit if the energy is high enough
        if (CumulativeEnergyDeposit > m_thresholdEnergyDeposit) {
	  DepositLowPtHit();
        }
        else { // energy is not enoug to ionize.
               // However, the track has ended and the energy is discarded and not added to the next step
               //G4cout << "reset due to end of track, discarding " << CumulativeEnergyDeposit/eV << " eV" << std::endl;
          ResetCumulativeVariables();
        }
      }
    }
    
    PreviousPostStepPosition = h.postPos();
    
    return true;
  }

  return true;
}

void TimeProjectionChamberSensitiveDetector::EndOfEvent(G4HCofThisEvent* HCE){
  // There might be one low Pt hit which has not been added to the collection.
  if (CumulativeEnergyDeposit > m_thresholdEnergyDeposit) {
    DepositLowPtHit(); 
  }  
}

void TimeProjectionChamberSensitiveDetector::DepositLowPtHit(){
  dd4hep::sim::Geant4TrackerHit* hit = new dd4hep::sim::Geant4TrackerHit(CurrentTrackID,CurrentPDGEncoding,CumulativeEnergyDeposit,CurrentGlobalTime);
  hit->cellID        = CurrentCopyNumber;
  hit->energyDeposit = CumulativeEnergyDeposit; // FIXME: also use the dedx
  hit->position      = CumulativeMeanPosition;
  hit->momentum      = CumulativeMeanMomentum;
  hit->length        = CumulativePathLength;
  m_lowPtHC->insert(hit);
  // reset the cumulative variables after positioning the hit
  ResetCumulativeVariables();
}

void TimeProjectionChamberSensitiveDetector::ResetCumulativeVariables(){
  CumulativeMeanPosition.SetXYZ(0.,0.,0.);// = dd4hep::Position(0.0,0.0,0.0);
  CumulativeMeanMomentum.SetXYZ(0.,0.,0.);// = dd4hep::sim::Momentum(0.0,0.0,0.0);
  CumulativeNumSteps = 0;
  CumulativeEnergyDeposit = 0;
  CumulativePathLength = 0;
}

void TimeProjectionChamberSensitiveDetector::CumulateLowPtStep(dd4hep::sim::Geant4StepHandler& h){
  //dd4hep::Position prePos    = h.prePos();
  //dd4hep::Position postPos   = h.postPos();
  //dd4hep::Position direction = postPos - prePos;
  dd4hep::Position meanPosition  = h.avgPosition();//0.5*(prePos+postPos);
  //double           hit_len   = direction.R();
  dd4hep::Position meanMomentum = 0.5*(h.preMom() + h.postMom());
  
  ++CumulativeNumSteps;    
  CumulativeMeanPosition   = ( (CumulativeMeanPosition*(CumulativeNumSteps-1)) + meanPosition ) / CumulativeNumSteps;
  CumulativeMeanMomentum   = ( (CumulativeMeanMomentum*(CumulativeNumSteps-1)) + meanMomentum ) / CumulativeNumSteps;
  CumulativeEnergyDeposit += h.deposit();
  CumulativePathLength    += h.stepLength();
  CurrentPDGEncoding       = h.trkPdgID();
  CurrentTrackID           = h.trkID();
  CurrentGlobalTime        = h.trkTime();
  CurrentCopyNumber        = h.preTouchable()->GetCopyNumber()/2+1;//getCellID( h.step );
  
  // needed for delta electrons: all particles causing hits have to be saved in the LCIO file 
  if ( CumulativeEnergyDeposit > m_thresholdEnergyDeposit ) {
    // This low Pt hit will eventually be saved, so set the flag to store the particle
    
    // writing the MC Particles can be turned on or off for the lowPt particles
    if(m_writeMCParticleForLowPtHits){
      // TODO: here set write flag for this track 

    }
  }
}
