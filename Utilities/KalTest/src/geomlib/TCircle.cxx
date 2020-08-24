//*************************************************************************
//* ====================
//*  TCircle Class
//* ====================
//*
//* (Description)
//*   A class to implement a circle object.
//* (Requires)
//*     TVCurve
//* (Provides)
//*     class TCircle
//* (Update Recored)
//*   2003/10/03  K.Fujii       Original version.
//*
//*************************************************************************
//
#include "TCircle.h"
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,15,0)
#include "TMath.h"
#endif
//_____________________________________________________________________
//  -----------------------------------
//  Circle Class
//  -----------------------------------

ClassImp(TCircle)

TCircle::TCircle(Double_t r, Double_t xc, Double_t yc)
       : fR(r), fXc(xc,yc)
{
}

Int_t TCircle::CalcXingPointWith(const TCircle  &c,
                                       TVector2  xx[],
                                       Double_t  eps) const
{
   TVector2 x12   = c.fXc - fXc;
   Double_t a     = x12.Mod2();
   Double_t sqa   = TMath::Sqrt(a);
   Double_t radd  = fR + c.fR;
   Double_t rsub  = fR - c.fR;
   Double_t arsub = TMath::Abs(rsub);
   if (sqa > radd+eps || sqa <= eps || sqa < arsub-eps) {
      // no intersection
      return 0;
   } else if (sqa > radd-eps) {
      // single intersection
      xx[0] = fXc + (fR/sqa)*x12;
      return 1;
   } else {
      // two intersection
      Double_t d = radd*rsub + a;
      a   = 1./sqa;
      d  *= 0.5*a;
      Double_t dp  = (fR+d)*(fR-d);
      if (dp <= 0.) dp = 0.;
      else          dp = TMath::Sqrt(dp)*a;
      d  *= a;
      TVector2 d1 = fXc + d*x12;
      TVector2 d2(dp*x12.Y(), -dp*x12.X());
      xx[0] = d1 + d2;
      xx[1] = d1 - d2;

      return 2;
   }
}

