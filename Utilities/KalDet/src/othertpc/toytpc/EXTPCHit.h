#ifndef EXTPCHIT_H
#define EXTPCHIT_H
//*************************************************************************
//* ================
//*  EXTPCHit Class
//* ================
//*
//* (Description)
//*   User defined hit class
//*   provides coordinate vector as defined by the MeasLayer
//* (Requires)
//*     TVTrackHit
//* (Provides)
//*     class EXTPCHit
//* (Update Recored)
//*   2009/11/23  K.Ikematsu   Derived from KalTest/examples/kaltest/
//*                                         hybrid/tpc/EXTPCHit.h
//*   2009/11/23  K.Ikematsu   Added GetSortingPolicy() method
//*   2009/11/23  K.Ikematsu   Added data member *fHitPtr
//*                            as a pointer to raw Hit object
//*
//* $Id: EXTPCHit.h,v 1.1 2010-03-11 15:07:01 fujiik Exp $
//*************************************************************************
//
#include "kaltest/KalTrackDim.h"
#include "kaltest/TVTrackHit.h"
#include "EXTPCMeasLayer.h"

class EXTPCHit : public TVTrackHit {

public:
  EXTPCHit(Int_t m = kMdim);

  EXTPCHit(const EXTPCMeasLayer &ms,
                 Double_t       *x,
                 Double_t       *dx,
                 Int_t           side,
                 Double_t        v,
           const TVector3       &xx,
                 Double_t        b,
                 Int_t           m = kMdim);

  EXTPCHit(const EXTPCMeasLayer &ms,
                 Double_t       *x,
                 Double_t       *dx,
                 Int_t           side,
                 Double_t        v,
           const void           *hitp,
                 Double_t        b,
                 Int_t           m = kMdim);

  virtual ~EXTPCHit();

  virtual TKalMatrix XvToMv(const TVector3 &xv, Double_t t0) const;
  virtual void       DebugPrint(Option_t *opt = "")          const;
  virtual Double_t   GetSortingPolicy()                      const;
  virtual Int_t      Compare(const TObject *obj)             const;
  virtual Bool_t     IsSortable()                            const;

  inline       Int_t     GetSide  () const { return fSide;   }
  inline       Double_t  GetVdrift() const { return fVdrift; }
  inline const void     *GetHitPtr() const { return fHitPtr; }
  inline const TVector3  GetExactX() const { return *fXXPtr; }

private:
  Int_t           fSide;    // (-1, +1) = (-z side, +z side)
  Double_t        fVdrift;  // drift veclocity
  const TVector3 *fXXPtr;   // pointer to exact hit
  const void     *fHitPtr;  // pointer to raw Hit object

  ClassDef(EXTPCHit, 1)  // EXTPC hit class
};
#endif
