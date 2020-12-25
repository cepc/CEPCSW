#ifndef TCYLINDER_H
#define TCYLINDER_H
//*************************************************************************
//* ====================
//*  TCylinder Class
//* ====================
//*
//* (Description)
//*   A class to implement a cylinder object.
//* (Requires)
//*     TVSurface
//* (Provides)
//*     class TCylinder
//* (Update Recored)
//*   2003/10/03  K.Fujii       Original version.
//*   2005/02/23  K.Fujii       Added GetSortingPolicy().
//*
//*************************************************************************
//
#include "TVSurface.h"
#include "TVector3.h"
#include "TMatrixD.h"

#include <cmath>

class TVTrack;

//_____________________________________________________________________
//  -----------------------------------
//  Cylinder Class
//  -----------------------------------

class TCylinder : public TVSurface {
public:
   TCylinder(Double_t r = 1., Double_t hlen  = 1., 
             Double_t xc = 0., Double_t yc = 0, Double_t zc = 0.)
           : fR(r), fHalfLen(hlen), fXc(xc,yc,zc) {}

   virtual ~TCylinder() {}

   virtual Int_t CalcXingPointWith(const TVTrack  &hel,
                                         TVector3 &xx,
                                         Double_t &phi,
                                         Int_t     mode,
                                         Double_t  eps  = 1.e-8) const;

   virtual Double_t CalcS   (const TVector3 &xx) const;
   virtual TMatrixD CalcDSDx(const TVector3 &xx) const;

   inline virtual       Bool_t     IsOnSurface(const TVector3 &xx) const;
   inline virtual       Bool_t     IsOutside  (const TVector3 &xx) const;

   inline virtual       Double_t   GetSortingPolicy()              const;

   inline virtual       Double_t   GetR      () const { return fR;    } 
   inline virtual const TVector3 & GetXc     () const { return fXc;   } 
   inline virtual       Double_t   GetLength () const;
   inline virtual       Double_t   GetZmin   () const;
   inline virtual       Double_t   GetZmax   () const;

private:
   Double_t fR;           // radius
   Double_t fHalfLen;     // half length
   TVector3 fXc;          // center
 
   ClassDef(TCylinder,1)      // cylinder class
};
//=======================================================
// inline functions
//=======================================================

Double_t TCylinder::GetLength () const 
{ 
   return 2*fHalfLen;
} 

Double_t TCylinder::GetZmin() const
{
   return fXc.Z() - fHalfLen;
}

Double_t TCylinder::GetZmax() const
{
   return fXc.Z() + fHalfLen;
}

Bool_t TCylinder::IsOnSurface(const TVector3 &xx) const
{
   return (xx.Z() >= GetZmin() && xx.Z() <= GetZmax()) && std::fabs( (xx-fXc).Perp() - fR ) < 1.e-6;
} 

Bool_t TCylinder::IsOutside(const TVector3 &xx) const
{
   Double_t r = (xx-fXc).Perp();
   Double_t z = xx.Z();
   return (r > fR || z < GetZmin() || z > GetZmax());
} 

Double_t TCylinder::GetSortingPolicy() const
{
   return GetR();
}
#endif

