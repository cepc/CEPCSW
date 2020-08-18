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
//*   2003/10/03  K.Fujii       Original version.  Currently fXc is 
//*                             supposed to be at the origin
//*   2009/05/30  K.Fujii       Now allow nonzero fXc.
//*
//*************************************************************************
//
#include <iostream>
#include "TCircle.h"
#include "TCylinder.h"
#include "TVTrack.h"

using namespace std;

//_____________________________________________________________________
//  -----------------------------------
//  Cylinder Class
//  -----------------------------------

ClassImp(TCylinder)

//_____________________________________________________________________
//  -----------------------------------
//  Calculate S
//  -----------------------------------
//
Double_t TCylinder::CalcS(const TVector3 &xx) const
{
   TVector3 xxc = xx - fXc;
   Double_t s   = xxc.X() * xxc.X() + xxc.Y() * xxc.Y() - fR * fR;
   return s;
}

//_____________________________________________________________________
//  -----------------------------------
//  Calculate (@S/@x)
//  -----------------------------------
//
TMatrixD TCylinder::CalcDSDx(const TVector3 &xx) const
{
   TVector3 xxc = xx - fXc;
   TMatrixD dsdx(1,3);
   dsdx(0,0) = 2.*xxc.X();
   dsdx(0,1) = 2.*xxc.Y();
   dsdx(0,2) = 0.;
   return dsdx;
}


//_____________________________________________________________________
//  -----------------------------------
//  Calculate crossing point with track
//  -----------------------------------
//
Int_t TCylinder::CalcXingPointWith(const TVTrack  &hel,
                                         TVector3 &xx,
                                         Double_t &phi,
                                         Int_t     mode,
                                         Double_t  eps) const
{
   // This assumes nonzero B field.
   //
   // Copy helix parameters to local variables.
   //
                                                                                
   Double_t dr  = hel.GetDrho();
   Double_t fi0 = hel.GetPhi0();
   Double_t cpa = hel.GetKappa();
   Double_t dz  = hel.GetDz();
   Double_t tnl = hel.GetTanLambda();
   TVector3 X0  = hel.GetPivot();
                                                                                
   //
   // Check if charge is nonzero.
   //
                                                                                
   Int_t    chg = (Int_t)TMath::Sign(1.1,cpa);
   if (!chg) {
      cerr << ">>>> Error >>>> TCylinder::CalcXingPointWith" << endl
           << "      Kappa = 0 is invalid for a helix "          << endl;
      return -1;
   }
                                                                                
   //
   // Project everything to XY plane and calculate crossing points.
   //
                                                                                
   Double_t rho  = hel.GetRho();
   Double_t rdr  = rho + dr;
   Double_t zdz  = X0.Z() + dz;
   Double_t csf0 = TMath::Cos(fi0);
   Double_t snf0 = TMath::Sin(fi0);
   Double_t xc  = X0.X() + rdr*csf0;
   Double_t yc  = X0.Y() + rdr*snf0;
   Double_t zc  = zdz;
                                                                                
   Double_t r   = TMath::Abs(rho);
   TCircle  c(r,xc,yc);
                                                                                
   Double_t rv  = GetR();
   Double_t xcv = GetXc().X();
   Double_t ycv = GetXc().Y();
   TCircle  cv(rv,xcv,ycv);
                                                                                
   TVector2 xxp[2];
   Int_t nx = c.CalcXingPointWith(cv, xxp);
   if (nx != 2) return 0;
                                                                                
   //
   // Crossing detected.
   //
                                                                                
   static const Double_t kPi     = TMath::Pi();
   static const Double_t kHalfPi = 0.5*TMath::Pi();
   static const Double_t kTwoPi  = 2.0*TMath::Pi();
                                                                                
   phi = 9999.;
   for (Int_t ix=0; ix<nx; ix++) {
      Double_t x   = xxp[ix].X() - xc;
      Double_t y   = xxp[ix].Y() - yc;
      Double_t dfi = TMath::ATan2(y,x) - fi0 - kHalfPi*(1+chg);
      if (!mode) {
         while (dfi < -kPi) dfi += kTwoPi;
         while (dfi >= kPi) dfi -= kTwoPi;
      } else {
         Int_t sign = (mode > 0 ? +1 : -1); // (+1,-1) = (fwd,bwd)
         while (dfi <  0.)     dfi += kTwoPi;
         while (dfi >= kTwoPi) dfi -= kTwoPi;
         if (sign*chg > 0) dfi -= kTwoPi;
      }
      if (TMath::Abs(dfi) < TMath::Abs(phi)) {
         phi = dfi;
         xx.SetXYZ(xxp[ix].X(), xxp[ix].Y(), zc - rho*tnl*phi);
      }
   }
   return (IsOnSurface(xx) ? 1 : 0);
}
