#ifndef TVTRACK_H
#define TVTRACK_H
//*************************************************************************
//* ===============
//*  TVTrack Class
//* ===============
//*
//* (Description)
//*   This is the base class for various tracks.
//* (Requires)
//*     TVCurve;
//* (Provides)
//*     class TVTrack
//* (Update Recored)
//*   2003/10/24  K.Fujii       Original version.
//*   2005/08/14  K.Fujii       Added IsInB().
//*
//*************************************************************************
//

#include "TVector3.h"
#include "TMatrixD.h"

#if 1
#include "TMath.h"
#include "TCollection.h"
#endif

#include "TVCurve.h"

//_____________________________________________________________________
//  -----------------------------------
//  Base Class for any track
//  -----------------------------------

class TVTrack : public TVCurve {
public:
   // Ctors and Dtor

   TVTrack(Double_t dr    = 0.,
           Double_t phi0  = 0.,
           Double_t kappa = 1.e-5,
           Double_t dz    = 0.,
           Double_t tanl  = 0.,
           Double_t x0    = 0.,
           Double_t y0    = 0.,
           Double_t z0    = 0.,
           Double_t b     = 30.);

   TVTrack(const TMatrixD &a, const TVector3 & x0, Double_t b = 30.);

   virtual ~TVTrack() {}

   // Utility methods

   virtual void MoveTo(const TVector3 &x0to,
                             Double_t &fid,
                             TMatrixD *F = 0,
                             TMatrixD *C = 0) = 0;

   virtual void MoveTo(const TVector3 &x0to,
                             Double_t &fid,
                             TMatrixD &C)
   {
      Int_t sdim = C.GetNrows();
      TMatrixD F(sdim, sdim);
      MoveTo(x0to,fid,&F,&C);
   }

   inline virtual void ScatterBy(Double_t dphi, Double_t dtnl)
   {
      fPhi0  += dphi;
      fKappa += (fKappa*fTanL/(1+fTanL*fTanL)) * dtnl;
      fTanL  += dtnl;
   }

   inline virtual void  PutInto(TMatrixD &sv) const
   {
      sv(0,0) = fDrho;
      sv(1,0) = fPhi0;
      sv(2,0) = fKappa;
      sv(3,0) = fDz;
      sv(4,0) = fTanL;
   }

   virtual TVector3 CalcXAt   (Double_t phi) const = 0;
   virtual TMatrixD CalcDxDa  (Double_t phi) const = 0;
   virtual TMatrixD CalcDxDphi(Double_t phi) const = 0;
   virtual void     CalcDapDa (Double_t fid,
                               Double_t dr,
                               Double_t drp,
                               TMatrixD &F)  const = 0;

   // Getters

   inline virtual       Double_t   GetDrho     () const { return fDrho;  }
   inline virtual       Double_t   GetPhi0     () const { return fPhi0;  }
   inline virtual       Double_t   GetKappa    () const { return fKappa; }
   inline virtual       Double_t   GetDz       () const { return fDz;    }
   inline virtual       Double_t   GetTanLambda() const { return fTanL;  }
   inline virtual const TVector3 & GetPivot    () const { return fX0;    }
   inline virtual       Double_t   GetRho      () const { return fAlpha/fKappa; }
   inline virtual       Double_t   GetPtoR     () const { return fAlpha; }

   // Setters

   inline virtual void  SetTo(const TMatrixD &sv, const TVector3 &x0)
   {
      fDrho  = sv(0,0);
      fPhi0  = sv(1,0);
      fKappa = sv(2,0);
      fDz    = sv(3,0);
      fTanL  = sv(4,0);
      fX0    = x0;
   }

   inline virtual void SetMagField(Double_t b)
   {
     // // units: mm, sec, Tesla
     if (b != 0.) fAlpha = kGiga/kLightVelocity*1000./b;
     // units: cm, sec, kGaus
     //if (b != 0.) fAlpha = kGiga/kLightVelocity*100./(b/10);
     else         fAlpha = kInfinity;
   }
	
   inline virtual Bool_t IsInB() const { return fAlpha < kInfinity ? kTRUE : kFALSE; }

protected:
   Double_t fDrho;      // drho
   Double_t fPhi0;      // phi0
   Double_t fKappa;     // kappa
   Double_t fDz;        // dz
   Double_t fTanL;      // tanl
   TVector3 fX0;        // pivot
   Double_t fAlpha;     // alpha

#if __GNUC__ < 4 && !defined(__STRICT_ANSI__)
   static const Double_t kLightVelocity = 2.99792458e8; //! light velocity [m/sec]
   static const Double_t kGiga          = 1.0e9;        //! Giga = 10^{9}
   static const Double_t kInfinity      = 1.e+20;       //! infinity
#else
   static const Double_t kLightVelocity;                //! light velocity [m/sec]
   static const Double_t kGiga;                         //! Giga = 10^{9}
   static const Double_t kInfinity;                     //! infinity
#endif

   ClassDef(TVTrack,1)      // Base class for any track
};

#endif
