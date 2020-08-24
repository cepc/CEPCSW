#ifndef THELICALTRACK_H
#define THELICALTRACK_H
//*************************************************************************
//* ====================
//*  THelicalTrack Class
//* ====================
//*
//* (Description)
//*   A class to implement a helical track object.
//*
//*   A helix is parametrized in the standard way:
//*
//*      x = x0 + drho * cos(phi0) + rho * (cos(phi0) - cos(phi0 + phi))
//*      y = y0 + drho * sin(phi0) + rho * (sin(phi0) - cos(sin0 + phi))
//*      z = z0 + dz               - rho * tan(lambda) * phi
//*
//*   with 
//*
//*      rho = alpha/kappa
//*
//* (Requires)
//*     TVTrack, TCylinder, TCircle, TVector3, TMatrixD
//* (Provides)
//*     class THelicalTrack
//* (Update Recored)
//*   2003/10/03  K.Fujii       Original version.
//*
//*************************************************************************
//

#include "TVTrack.h"

//_____________________________________________________________________
//  -----------------------------------
//  Helical Track Class
//  -----------------------------------

class THelicalTrack : public TVTrack {
public:

   // Ctors and Dtor

   THelicalTrack(Double_t dr    = 0.,
                 Double_t phi0  = 0.,
                 Double_t kappa = 1.e-5,
                 Double_t dz    = 0.,
                 Double_t tanl  = 0.,
                 Double_t x0    = 0.,
                 Double_t y0    = 0.,
                 Double_t z0    = 0.,
                 Double_t b     = 30.);

   THelicalTrack(const TMatrixD &a, const TVector3 &x0, Double_t b = 30.);
   THelicalTrack(const TVector3 &x1, const TVector3 &x2, const TVector3 &x3,
                 Double_t b = 30., Bool_t dir = kIterForward);

   virtual ~THelicalTrack() {}

   // Utility methods

   virtual void MoveTo(const TVector3 &x0to,
                             Double_t &fid,
                             TMatrixD *F = 0,
                             TMatrixD *C = 0);

   TVector3 CalcXAt   (Double_t phi) const;
   TMatrixD CalcDxDa  (Double_t phi) const;
   TMatrixD CalcDxDphi(Double_t phi) const;
   void     CalcDapDa (Double_t fid,
                       Double_t dr,
                       Double_t drp,
                       TMatrixD &F)  const;

private:
   void CalcStartHelix(const TVector3 &x1,
                       const TVector3 &x2,
                       const TVector3 &x3,
                             Bool_t    dir = kIterForward);

private:

   ClassDef(THelicalTrack,1)      // circle class
};

#endif

