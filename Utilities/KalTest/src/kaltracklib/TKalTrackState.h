#ifndef TKALTRACKSTATE_H
#define TKALTRACKSTATE_H
//*************************************************************************
//* ======================
//*  TKalTrackState Class
//* ======================
//*
//* (Description)
//*   Track state vector class used in Kalman Filter.
//* (Requires)
//*     TVKalState
//* (Provides)
//*     class TKalTrackState
//* (Update Recored)
//*   2003/09/30  Y.Nakashima       Original version.
//*   2005/02/23  A.Yamaguchi       Added CalcDapDa method.
//*   2005/02/XX  A.Yamaguchi       Moved CalcDapDa method to THelicalTrack.
//*   2005/08/13  K.Fujii           Removed CalcProcessNoise method.
//*   2010/04/06  K.Fujii           Modified MoveTo to allow a 1-dim hit,
//*                                 for which pivot is at the xpected hit.
//*
//*************************************************************************

#include "TVKalState.h"         // from KalLib
#include "THelicalTrack.h"      // from GeomLib
#include "TStraightTrack.h"     // from GeomLib
#include "KalTrackDim.h"        // from KalTrackLib

class TKalTrackSite;

//_________________________________________________________________________
//  -----------------------------------
//  Base Class for Kalman state vector
//  -----------------------------------
//

class TKalTrackState : public TVKalState {
                                                                                
public:
                                                                                
   // Ctors and Dtor
                                                                                
   TKalTrackState(Int_t p = kSdim);
   TKalTrackState(const TKalMatrix &sv, Int_t type = 0, Int_t p = kSdim);
   TKalTrackState(const TKalMatrix &sv, const TKalMatrix &c, 
                        Int_t type = 0, Int_t p = kSdim);
   TKalTrackState(const TKalMatrix &sv, const TVKalSite &site, 
                        Int_t type = 0, Int_t p = kSdim);
   TKalTrackState(const TKalMatrix &sv, const TKalMatrix &c,
                  const TVKalSite &site, Int_t type = 0, Int_t p = kSdim);
   virtual ~TKalTrackState() {}
                                                                                
   // Implementation of paraent class pure virtuals
                                                                                
   TKalTrackState * MoveTo(TVKalSite  &to, 
                           TKalMatrix &F, 
                           TKalMatrix *QPtr = 0) const;
   TKalTrackState & MoveTo(TVKalSite  &to, 
                           TKalMatrix &F, 
                           TKalMatrix &Q) const;
   void         DebugPrint() const;

   // Derived class methods

   THelicalTrack   GetHelix() const;
   TStraightTrack  GetLine () const;
   TVTrack        &CreateTrack() const;

private:

   TVector3 fX0;		// pivot

   ClassDef(TKalTrackState,1)      // sample state vector class
};
                                                                                
#endif
