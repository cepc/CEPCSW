#ifndef TCUTCONE_H
#define TCUTCONE_H
//*************************************************************************
//* ====================
//*  TCutCone Class
//* ====================
//*
//* (Description)
//*   A class to implement a conical surface object.
//* (Requires)
//*     TVSurface
//* (Provides)
//*     class TCutCone
//* (Update Recored)
//*   2012/01/19  K.Fujii       Original version derived from THYpe.
//*                             This class is implemented as an extreme
//*                             (fR0=0) case of THype, thereby containing
//*                             both -ve and +ve sides. If you want to
//*                             restrict them to one side, override
//*                             IsOnSurface(), etc.
//*************************************************************************
//
#include "TVSurface.h"
#include "TVector3.h"
#include "TMatrixD.h"

class TVTrack;

//_____________________________________________________________________
//  -----------------------------------
//  Hype Class
//  -----------------------------------

class TCutCone : public TVSurface {
public:
   TCutCone(Double_t z1 = 0.5,
            Double_t hlen = 1.0,
            Double_t tana = 0.1,
            Double_t xc = 0., 
            Double_t yc = 0, 
            Double_t zc = 0.)
      : fZ1(z1), 
        fHalfLen(hlen), 
        fXc(xc,yc,zc), 
        fTanA(tana)
   {
   }

   virtual ~TCutCone() {}

   virtual Double_t CalcS   (const TVector3 &xx) const;
   virtual TMatrixD CalcDSDx(const TVector3 &xx) const;

   inline virtual       Bool_t     IsOnSurface(const TVector3 &xx) const;
   inline virtual       Bool_t     IsOutside  (const TVector3 &xx) const;

   inline virtual       Double_t   GetSortingPolicy()              const;

   inline virtual       Double_t   GetZ1     () const { return fZ1;   } 
   inline virtual const TVector3 & GetXc     () const { return fXc;   } 
   inline virtual       Double_t   GetTanA   () const { return fTanA; }
   inline virtual       Double_t   GetLength () const;
   inline virtual       Double_t   GetZmin   () const;
   inline virtual       Double_t   GetZmax   () const;

private:
   Double_t fZ1;          // z position of the front face
   Double_t fHalfLen;     // half length (cone length from the apex)
   TVector3 fXc;          // center
   Double_t fTanA;        // tan(half cone angle)
#if __GNUC__ < 4 && !defined(__STRICT_ANSI__)
   static const Double_t kTol = 1.e-5; // tolerance
#else
   static const Double_t kTol; // tolerance
#endif
 
   ClassDef(TCutCone,1)      // hype class
};
//=======================================================
// inline functions
//=======================================================

Double_t TCutCone::GetLength () const 
{ 
   return 2*fHalfLen;
} 

Double_t TCutCone::GetZmin() const
{
   return TMath::Min(fXc.Z() - fHalfLen, fXc.Z() + fHalfLen);
}

Double_t TCutCone::GetZmax() const
{
   return TMath::Max(fXc.Z() - fHalfLen, fXc.Z() + fHalfLen);
}

Bool_t TCutCone::IsOnSurface(const TVector3 &xx) const
{
   TVector3 xxc = xx - fXc;
   Double_t r   = xxc.Perp();
   Double_t z   = xxc.Z();
   Double_t s   = r*r - fTanA*fTanA*z*z;
   
   return (TMath::Abs(s) < kTol && xx.Z() >= GetZmin() && xx.Z() <= GetZmax());
} 

Bool_t TCutCone::IsOutside(const TVector3 &xx) const
{
   Double_t r  = (xx-fXc).Perp();
   Double_t z  = xx.Z();
   Double_t R2 =fTanA*fTanA*z*z;
   
   return (r*r > R2 || z < GetZmin() || z > GetZmax());
} 

Double_t TCutCone::GetSortingPolicy() const
{
   return GetZ1()*GetTanA();
}
#endif

