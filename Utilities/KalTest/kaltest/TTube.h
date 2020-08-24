#ifndef TTUBE_H
#define TTUBE_H
//*************************************************************************
//* ====================
//*  TTube Class
//* ====================
//*
//* (Description)
//*   A class to implement a tube object.
//* (Requires)
//*     TVSolid
//* (Provides)
//*     class TTube
//* (Update Recored)
//*   2003/10/03  K.Fujii       Original version.
//*
//*************************************************************************
//
#include "TVSolid.h"
#include "TVector3.h"

class TVTrack;

//_____________________________________________________________________
//  -----------------------------------
//  Cylinder Class
//  -----------------------------------

class TTube : public TVSolid {
public:
   TTube(Double_t rin =  0., Double_t rout = 1., Double_t hlen  = 1., 
             Double_t xc = 0., Double_t yc = 0, Double_t zc = 0.)
           : fRin(rin), fRout(rout), fHalfLen(hlen), fXc(xc,yc,zc) {}

   virtual ~TTube() {}

   virtual Int_t CalcXingPointWith(const TVTrack  &hel,
                                         Double_t &phi,
                                         TVector3 &xx,
                                         Int_t     mode,
                                         Double_t  eps = 1.e-8)   const;

   inline virtual       Bool_t     IsOnBarrel(const TVector3 &xx) const;
   inline virtual       Bool_t     IsOutside (const TVector3 &xx) const;

   inline virtual       Double_t   GetRin    () const { return fRin;  } 
   inline virtual       Double_t   GetRout   () const { return fRout; } 
   inline virtual const TVector3 & GetXc     () const { return fXc;   } 
   inline virtual       Double_t   GetLength () const;
   inline virtual       Double_t   GetZmin   () const;
   inline virtual       Double_t   GetZmax   () const;

private:
   Double_t fRin;         // inner radius
   Double_t fRout;        // outer radius
   Double_t fHalfLen;     // half length
   TVector3 fXc;          // center
 
   ClassDef(TTube,1)      // TTube class
};
//=======================================================
// inline functions
//=======================================================

Double_t TTube::GetLength () const 
{ 
   return 2*fHalfLen;
} 

Double_t TTube::GetZmin() const
{
   return fXc.Z() - fHalfLen;
}

Double_t TTube::GetZmax() const
{
   return fXc.Z() + fHalfLen;
}

Bool_t TTube::IsOnBarrel(const TVector3 &xx) const
{
   return (xx.Z() >= GetZmin() && xx.Z() <= GetZmax());
} 

Bool_t TTube::IsOutside(const TVector3 &xx) const
{
   Double_t r = (xx-fXc).Perp();
   Double_t z = xx.Z();
   return (r < fRin || r > fRout || z < GetZmin() || z > GetZmax());
} 
#endif

