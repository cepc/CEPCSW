#ifndef TVKALSTATE_H
#define TVKALSTATE_H
//*************************************************************************
//* ====================
//*  TVKalState Class
//* ====================
//*
//* (Description)
//*   This is the base class for a state vector used in Kalman Filter.
//* (Requires)
//* 	TKalMatrix
//* (Provides)
//* 	class TVKalState
//* (Update Recored)
//*   2003/09/30  K.Fujii	Original version.
//*
//*************************************************************************
//
#include "TKalMatrix.h"
class TVKalSite;
//_____________________________________________________________________
//  -----------------------------------
//  Base Class for Kalman state vector
//  -----------------------------------
//
class TVKalState : public TKalMatrix { 
public:

   // Ctors and Dtor

   TVKalState(Int_t type = 0, Int_t p = 6);
   TVKalState(const TKalMatrix &sv, Int_t type = 0, Int_t p = 6);
   TVKalState(const TKalMatrix &sv, const TKalMatrix &c, 
              Int_t type = 0, Int_t p = 6);
   TVKalState(const TKalMatrix &sv, const TVKalSite &site, 
              Int_t type = 0, Int_t p = 6);
   TVKalState(const TKalMatrix &sv, const TKalMatrix &c,
              const TVKalSite &site, Int_t type = 0, Int_t p = 6);
   virtual ~TVKalState() {}

   // Pure virtuals to be implemented in derived classes
   //
   // MoveTo should calculate
   //    a:  predicted state vector      : a^k-1_k = f_k-1(a_k-1)
   //    F:  propagator derivative       : F_k-1   = (@f_k-1/@a_k-1)
   //    Q:  process noise from k-1 to k : Q_k-1)
   // and return a^k-1_k.
   //

   virtual TVKalState * MoveTo(TVKalSite  &to,
                               TKalMatrix &F,
                               TKalMatrix *QPtr = 0) const = 0;
   virtual TVKalState & MoveTo(TVKalSite  &to,
                               TKalMatrix &F,
                               TKalMatrix &Q) const = 0;

   virtual void         DebugPrint() const = 0;

   virtual void         Propagate(TVKalSite &to); // calculates f, F, and Q
 
   
   // Getters

   inline virtual Int_t GetDimension                () const { return GetNrows(); }
   inline virtual const TVKalSite  & GetSite        () const { return *fSitePtr; }
   inline virtual const TKalMatrix & GetCovMat      () const { return fC; }
   inline virtual const TKalMatrix & GetProcNoiseMat() const { return fQ; }
   inline virtual const TKalMatrix & GetPropMat     (const Char_t *t = "") const { return (t[0] == 'T' ? fFt : fF); } 

   // Setters

   inline virtual void SetStateVec    (const TKalMatrix &c) { TMatrixD::operator=(c); }
   inline virtual void SetCovMat      (const TKalMatrix &c) { fC       = c; }
   inline virtual void SetProcNoiseMat(const TKalMatrix &q) { fQ       = q; }
   inline virtual void SetSitePtr     (TVKalSite  *s)       { fSitePtr = s; }

private:
   
   // private data members -------------------------------------------

   Int_t       fType;    // (0,1,2,3) = (uninited,predicted,filtered,smoothed)
   TVKalSite  *fSitePtr; // pointer to corresponding KalSite
   TKalMatrix  fF;       // propagator matrix to next site (F = @f/@a)
   TKalMatrix  fFt;      // transposed propagator matrix (F^T = (@f/@a)^T)
   TKalMatrix  fQ;       // process noise from this to the next sites
   TKalMatrix  fC;       // covariance matrix
  
   ClassDef(TVKalState,1)      // Base class for state vector objects
};
#endif
