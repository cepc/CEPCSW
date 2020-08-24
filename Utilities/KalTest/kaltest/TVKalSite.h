#ifndef TVKALSITE_H
#define TVKALSITE_H
//*************************************************************************
//* ===================
//*  TVKalSite Class
//* ===================
//*
//* (Description)
//*   This is the base class for measurement vector used by Kalman filter.
//* (Requires)
//* 	TKalMatrix
//* (Provides)
//* 	class TVKalSite
//* (Update Recored)
//*   2003/09/30  K.Fujii	Original version.
//*   2005/02/23  A.Yamaguchi	Added getter and setter for a new static
//*                             data member, fgKalSysPtr.
//*   2005/08/25  A.Yamaguchi	Removed getter and setter for a new static
//*                             data member, fgKalSysPtr.
//*   2009/06/18  K.Fujii       Implement inverse Kalman filter
//*
//*************************************************************************
//
#include "TObjArray.h"

#include "TAttLockable.h"
#include "TKalMatrix.h"
#include "TVKalState.h"

class TVKalSystem;

//_____________________________________________________________________
//  ------------------------------
//  Base Class for Kalman measurement vector
//  ------------------------------
//
class TVKalSite : public TObjArray, public TAttLockable {
friend class TVKalSystem;

public:
   enum EStType { kPredicted = 0,
                  kFiltered,
                  kSmoothed,
                  kInvFiltered };
public:
   // Ctors and Dtor

   TVKalSite(Int_t m = 2, Int_t p = 6);
   virtual ~TVKalSite();
             
   // Utility Methods

   virtual Int_t   CalcExpectedMeasVec  (const TVKalState &a,
                                               TKalMatrix &m) = 0;
   virtual Int_t   CalcMeasVecDerivative(const TVKalState &a,
                                               TKalMatrix &H) = 0;
   virtual Bool_t  IsAccepted() = 0;

   virtual void    DebugPrint() const = 0;

   virtual Bool_t  Filter();

   virtual void    Smooth(TVKalSite &pre);

   virtual void    InvFilter();

   inline  void    Add(TObject *obj);
   
   // Getters

   inline virtual Int_t        GetDimension() const { return fM.GetNrows(); }
   inline virtual TVKalState & GetCurState ()       { return *fCurStatePtr; }
   inline virtual TVKalState & GetCurState () const { return *fCurStatePtr; }
   inline virtual TVKalState & GetState (EStType t);
   inline virtual TKalMatrix & GetMeasVec      ()   { return fM;            }
   inline virtual TKalMatrix & GetMeasNoiseMat ()   { return fV;            }
   inline virtual TKalMatrix & GetResVec       ()   { return fResVec;       }
   inline virtual TKalMatrix & GetCovMat       ()   { return fR;            }
   inline virtual Double_t     GetDeltaChi2() const { return fDeltaChi2;    }
          virtual TKalMatrix   GetResVec (EStType t);

private:
   // Private utility methods

   virtual TVKalState & CreateState(const TKalMatrix &sv, Int_t type = 0) = 0;
   virtual TVKalState & CreateState(const TKalMatrix &sv, const TKalMatrix &c,
                                    Int_t type = 0) = 0;

private:
   
   // private data member -------------------------------------------

   TVKalState    *fCurStatePtr; // pointer to current best state
   TKalMatrix     fM;           // measurement vector: M(m,1)
   TKalMatrix     fV;           // noise matrix: M(m,m)
   TKalMatrix     fH;           // H = (@h/@a): M(m,p)
   TKalMatrix     fHt;          // H^t = (@h/@a)^t: M(p,m)
   TKalMatrix     fResVec;      // m - h(a): M(m,1)
   TKalMatrix     fR;           // covariance matrix: M(m,m)
   Double_t       fDeltaChi2;   // chi2 increment

   ClassDef(TVKalSite,1)      // Base class for measurement vector objects
};

//=======================================================
// inline functions
//=======================================================

void TVKalSite::Add(TObject *obj)
{
   TObjArray::Add(obj);
   fCurStatePtr = static_cast<TVKalState *>(obj);
   fCurStatePtr->SetSitePtr(this);
}

TVKalState & TVKalSite::GetState(TVKalSite::EStType t)
{
   TVKalState *ap = 0;
   if (t >= 0 && t < GetEntries()) {
      ap = static_cast<TVKalState *>(UncheckedAt(t));
   }
   return *ap;
}
#endif
