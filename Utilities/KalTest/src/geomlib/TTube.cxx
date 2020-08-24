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
#include <iostream>
#include "TCircle.h"
#include "TTube.h"
#include "TVTrack.h"

using namespace std;

//_____________________________________________________________________
//  -----------------------------------
//  Tube Class
//  -----------------------------------

ClassImp(TTube)

//_____________________________________________________________________
//  -----------------------------------
//  Calculate crossing point with helix
//  -----------------------------------
//
Int_t TTube::CalcXingPointWith(const TVTrack  &hel,
                                     Double_t &phi,
                                     TVector3 &xx,
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
      cerr << ">>>> Error >>>> TTube::CalcXingPointWith" << endl
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
                                                                                
   Double_t rv  = GetRout();
   Double_t xcv = GetXc().X();
   Double_t ycv = GetXc().Y();
   TCircle  cv(rv,xcv,ycv);
                                                                                
   TVector2 xxp[2];
   Int_t nx = c.CalcXingPointWith(cv, xxp);
                                                                                
   //
   // Switch on the number of crossing points.
   //
                                                                                
   static const Double_t kPi     = TMath::Pi();
   static const Double_t kHalfPi = 0.5*TMath::Pi();
   static const Double_t kTwoPi  = 2.0*TMath::Pi();
                                                                                
   Int_t iret = 0;
   switch (nx) {
      case 2:  // track hitting barrel part
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
               xx.SetXYZ(xxp[ix].X(), xxp[ix].Y(), zc - rho*tnl*dfi);
            }
         }
         if (IsOnBarrel(xx)) return 1;
      default: // track hitting end cap part
         if (TMath::Abs(tnl) < 0.1) {
            return 0;
         } else if (tnl < 0.) {
            xx.SetZ(GetZmin());
            iret = 2;
         } else {
            xx.SetZ(GetZmax());
            iret = 3;
         }
         phi = (zdz-xx.Z())/rho/tnl;
         if (TMath::Abs(phi) > kTwoPi) {
            return 0;
         } else {
            xx.SetX(xc - rho*TMath::Cos(phi+fi0));
            xx.SetY(yc - rho*TMath::Sin(phi+fi0));
            if (IsOutside(xx)) return 0;
         }
         break;
   }
   return iret;
}
