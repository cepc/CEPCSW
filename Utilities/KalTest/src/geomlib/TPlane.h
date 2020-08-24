#ifndef TPLANE_H
#define TPLANE_H
//*************************************************************************
//* ====================
//*  TPlane Class
//* ====================
//*
//* (Description)
//*   A class to implement a flat plane object.
//* (Requires)
//*     TVSurface
//* (Provides)
//*     class TPlane
//* (Update Recored)
//*   2004/10/30  A.Yamaguchi       Original version.
//*   2005/02/23  K.Fujii           Added GetSortingPolicy().
//*
//*************************************************************************
//
#include "TVSurface.h"
#include "TVector3.h"
#include "TMatrixD.h"

class TVTrack;

//_____________________________________________________________________
//  -----------------------------------
//  Plane Class
//  -----------------------------------

class TPlane : public TVSurface {
public:
   TPlane();
   TPlane(const TVector3 &xc);
   TPlane(const TVector3 &xc, const TVector3 &n);

   virtual ~TPlane() {}

   virtual Double_t CalcS   (const TVector3 &xx) const;
   virtual TMatrixD CalcDSDx(const TVector3 &xx) const;

   inline virtual const TVector3 & GetXc     () const { return fXc;     } 
   inline virtual const TVector3 & GetNormal () const { return fNormal; } 
   inline virtual       Bool_t     IsOnSurface(const TVector3 &xx) const;
   inline virtual       Bool_t     IsOutside  (const TVector3 &xx) const;

   inline virtual       Double_t   GetSortingPolicy()              const;


private:
   TVector3 fXc;          // center
   TVector3 fNormal;      // normal
 
   ClassDef(TPlane,1)      // plane class
};
//=======================================================
// inline functions
//=======================================================


Bool_t TPlane::IsOnSurface(const TVector3 &xx) const
{
#if 0
   return (xx - fXc) * fNormal == 0. ? kTRUE : kFALSE; 
#else
   return kTRUE;
#endif
} 

Bool_t TPlane::IsOutside(const TVector3 &xx) const
{
   return (xx - fXc) * fNormal > 0. ? kTRUE : kFALSE; 
} 

Double_t TPlane::GetSortingPolicy() const
{
   return TMath::Abs(fXc*fNormal.Unit());
}
#endif

