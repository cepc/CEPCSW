//*************************************************************************
//* ====================
//*  THelicalTrack Class
//* ====================
//*
//* (Description)
//*   A class to implement a helical track object.
//* (Requires)
//*     TVTrack, TCylinder, TCircle, TVector3, TMatrixD
//* (Provides)
//*     class THelicalTrack
//* (Update Recored)
//*   2003/10/03  K.Fujii       Original version.
//*
//*************************************************************************
//

#include <iostream>
#include "THelicalTrack.h"

using namespace std;

//_____________________________________________________________________
//  -----------------------------------
//  Helical Track Class
//  -----------------------------------
//
//_____________________________________________________________________
//  --------------
//  Ctors and Dtor
//  --------------
//

THelicalTrack::THelicalTrack(Double_t dr,
                             Double_t phi0,
                             Double_t kappa,
                             Double_t dz,
                             Double_t tanl,
                             Double_t x0,
                             Double_t y0,
                             Double_t z0,
                             Double_t b)
             : TVTrack(dr,phi0,kappa,dz,tanl, x0,y0,z0, b)
{
}

THelicalTrack::THelicalTrack(const TMatrixD &a, const TVector3 &x0, Double_t b)
             : TVTrack(a, x0, b)
{
}

THelicalTrack::THelicalTrack(const TVector3 &x1, 
                             const TVector3 &x2,
                             const TVector3 &x3,
                                   Double_t  b,
                                   Bool_t    dir)
{
   SetMagField(b);
   CalcStartHelix(x1,x2,x3,dir);
}

//_____________________________________________________________________
//  ----------------
//  Utility methods
//  ----------------

void THelicalTrack::MoveTo(const TVector3 &xv0to, // new pivoit
                                 Double_t &fid,   // deflection angle
                                 TMatrixD *FPtr,  // propagator matrix
                                 TMatrixD *CPtr)  // covariance matrix
{
   // ---------------------------------------------------
   // (0) Preparation
   // ---------------------------------------------------
   //   Define some numerical constants.
   //

   static const Double_t kPi    = TMath::Pi();
   static const Double_t kTwoPi = 2.0*kPi;

   //   Copy helix parmeters to local variables

   Double_t dr    = fDrho;
   Double_t fi0   = fPhi0;
   while (fi0 < 0.)      fi0 += kTwoPi;
   while (fi0 > kTwoPi)  fi0 -= kTwoPi;
   Double_t cpa   = fKappa;
   Double_t dz    = fDz;
   Double_t tnl   = fTanL;

   Double_t x0    = fX0.X();
   Double_t y0    = fX0.Y();
   Double_t z0    = fX0.Z();
   Double_t xv    = xv0to.X();
   Double_t yv    = xv0to.Y();
   Double_t zv    = xv0to.Z();

   // ---------------------------------------------------
   // (1) Calculate a' = f_k-1(a_k-1)
   // ---------------------------------------------------
   //        a' = (dr', fi0', cpa', dz', tnl')
   //        a  = (dr , fi0 , cpa , dz , tnl )
   //

   Double_t r     = fAlpha/cpa;
   Double_t rdr   = r + dr;
   Double_t csf0  = TMath::Cos(fi0);
   Double_t snf0  = TMath::Sqrt(TMath::Max(0.0, (1.0-csf0)*(1.0+csf0)));
   if (fi0 > kPi) snf0 = -snf0;

   Double_t xc    = x0 + rdr*csf0;
   Double_t yc    = y0 + rdr*snf0;
   Double_t fi0p  = 0.;

   if (cpa > 0.) fi0p = TMath::ATan2((yc-yv),(xc-xv));
   if (cpa < 0.) fi0p = TMath::ATan2((yv-yc),(xv-xc));
   while (fi0p < 0.)      fi0p += kTwoPi;
   while (fi0p > kTwoPi)  fi0p -= kTwoPi;

   Double_t csf   = TMath::Cos(fi0p);
   Double_t snf   = TMath::Sqrt(TMath::Max(0.0, (1.0-csf)*(1.0+csf)));
   if (fi0p > kPi) snf = -snf;

   Double_t anrm  = 1.0/TMath::Sqrt(csf*csf+snf*snf);
            csf  *= anrm;
            snf  *= anrm;
   Double_t csfd  = csf*csf0 + snf*snf0;
   Double_t snfd  = snf*csf0 - csf*snf0;

   fid   = fi0p - fi0;
   while (fid < 0)      fid += kTwoPi;
   while (fid > kTwoPi) fid -= kTwoPi;
   if    (fid > kPi)    fid -= kTwoPi;

   Double_t drp   = (xc-xv)*csf + (yc-yv)*snf - r;
   Double_t dzp   = z0 - zv + dz - r*tnl*fid;
  
  
#if 0 //======== this code should not be necessary with a bug fix in TKalDetCradle =====================

  // make sure that the helix really moves to the closest point to the reference point
  // use protective_counter to ensure we don't enter an infinate loop
  
  if (tnl > 1.0e-10) { // protect against unreasonably small values of tnl, as in the case special case of a circle this will lead to division by zero below 
    
    double phi_to_ref =  dzp/(-r*tnl);
    int protective_counter = 0;
    
    while ((phi_to_ref < -kPi) && protective_counter < 1000 ) {
      phi_to_ref += kTwoPi;
      ++protective_counter;
    }
    
    protective_counter = 0;
    
    while ((phi_to_ref >  kPi) && protective_counter < 1000 ) {
      phi_to_ref -= kTwoPi;
      ++protective_counter;
    }
    
    double phi_correction = dzp/(-r*tnl) - phi_to_ref;
    
    dzp += phi_correction*r*tnl;

  }
  
#endif // =============================================================================================


   TMatrixD av(5,1);
   av(0,0) = drp;
   av(1,0) = fi0p;
   av(2,0) = cpa;
   av(3,0) = dzp;
   av(4,0) = tnl;

   THelicalTrack helto(av,xv0to);
   helto.fAlpha = fAlpha;

   if (!FPtr && !CPtr) {
      *this = helto;
      return;
   }

   TMatrixD Fdummy(5,5);
   TMatrixD &F = FPtr ? *FPtr : Fdummy;

   // ---------------------------------------------------
   // (2) Calculate @a'/@a = @a'/a = F_k-1
   // ---------------------------------------------------
   //        a' = (dr', fi0', cpa', dz', tnl')
   //        a  = (dr , fi0 , cpa , dz , tnl )
   //

   Double_t rdrpr = 1.0/(r+drp);
   Double_t rcpar = r/cpa;

   // @drho'/@a
   F(0,0) = csfd;
   F(0,1) = rdr*snfd;
   F(0,2) = rcpar*(1.0-csfd);
   F(0,3) = 0;
   F(0,4) = 0;

   // @phi0'/@a
   F(1,0) = -rdrpr*snfd;
   F(1,1) =  rdr*rdrpr*csfd;
   F(1,2) =  rcpar*rdrpr*snfd;
   F(1,3) =  0;
   F(1,4) =  0;

   // @kappa'/@a
   F(2,0) = 0;
   F(2,1) = 0;
   F(2,2) = 1;
   F(2,3) = 0;
   F(2,4) = 0;

   // @dz'/@a
   F(3,0) =  r*rdrpr*tnl*snfd;
   F(3,1) =  r*tnl*(1.0-rdr*rdrpr*csfd);
   F(3,2) =  rcpar*tnl*(fid-r*rdrpr*snfd);
   F(3,3) =  1;
   F(3,4) = -r*fid;

   // @tanl'/@a
   F(4,0) = 0;
   F(4,1) = 0;
   F(4,2) = 0;
   F(4,3) = 0;
   F(4,4) = 1;

   if (!CPtr) {
      *this = helto;
      return;
   }

   // ---------------------------------------------------
   // (3) Calculate C' = C^k-1_k
   // ---------------------------------------------------

   TMatrixD &C  = *CPtr;
   TMatrixD  Ft = TMatrixD(TMatrixD::kTransposed, F);
   TMatrixD  Cp = F * C * Ft;
   C = Cp;

   *this = helto;
}

TVector3 THelicalTrack::CalcXAt(Double_t phi) const
{
   Double_t csf0 = TMath::Cos(fPhi0);
   Double_t snf0 = TMath::Sin(fPhi0);
   Double_t snfd = TMath::Sin(fPhi0 + phi);
   Double_t csfd = TMath::Cos(fPhi0 + phi);
   Double_t rho  = fAlpha/fKappa;

   Double_t x    = fX0.X() + fDrho * csf0 + rho * (csf0 - csfd);
   Double_t y    = fX0.Y() + fDrho * snf0 + rho * (snf0 - snfd);
   Double_t z    = fX0.Z() + fDz          - rho * fTanL * phi;

#if 0
   std::cerr << "THelicalTrack::CalcXAt" << std::endl;
   std::cerr << " phi0 = " << fPhi0 << std::endl;
   std::cerr << " phi  = " << phi   << std::endl;
   std::cerr << " Drho = " << fDrho << std::endl;
   std::cerr << " rho  = " << rho   << std::endl;
   std::cerr << " x0   = " << fX0.X()     << std::endl;
   std::cerr << " y0   = " << fX0.Y()     << std::endl;
   std::cerr << " z0   = " << fX0.Z()     << std::endl;
   std::cerr << " x    = " << x     << std::endl;
   std::cerr << " y    = " << y     << std::endl;
   std::cerr << " z    = " << z     << std::endl;
#endif

   return TVector3(x,y,z);
}

TMatrixD THelicalTrack::CalcDxDa(Double_t phi) const
{
   Double_t fi0   = fPhi0;
   Double_t r     = fAlpha/fKappa;
   Double_t rcpar = r/fKappa;

   Double_t snf0 = TMath::Sin(fi0);
   Double_t csf0 = TMath::Cos(fi0);
   Double_t snfd = TMath::Sin(fi0 + phi);
   Double_t csfd = TMath::Cos(fi0 + phi);

   TMatrixD dxda(3,5);
   // @x/@a
   dxda(0,0) =  csf0;
   dxda(0,1) = -fDrho * snf0 - r * (snf0 - snfd);
   dxda(0,2) = -rcpar * (csf0 - csfd);
   dxda(0,3) =  0;
   dxda(0,4) =  0;

   // @y/@a
   dxda(1,0) =  snf0;
   dxda(1,1) =  fDrho * csf0 + r *(csf0 - csfd);
   dxda(1,2) = -rcpar * (snf0 - snfd);
   dxda(1,3) =  0;
   dxda(1,4) =  0;

   // @z/@a
   dxda(2,0) = 0;
   dxda(2,1) = 0;
   dxda(2,2) = rcpar * phi * fTanL;
   dxda(2,3) = 1;
   dxda(2,4) = - r * phi;

   return dxda;
}

TMatrixD THelicalTrack::CalcDxDphi(Double_t phi) const
{
   Double_t r    = fAlpha/fKappa;

   Double_t snfd = TMath::Sin(fPhi0 + phi);
   Double_t csfd = TMath::Cos(fPhi0 + phi);

   TMatrixD dxdphi(3,1);
   dxdphi(0,0) =  r * snfd;
   dxdphi(1,0) = -r * csfd;
   dxdphi(2,0) = -r * fTanL;

   return dxdphi;
}

void THelicalTrack::CalcStartHelix(const TVector3 &x1,
                                   const TVector3 &x2,
                                   const TVector3 &x3,
                                         Bool_t    dir)
{
   static const TVector3 ez(0., 0., 1.);
   TVector3 x12 = x2 - x1;
   TVector3 x13 = x3 - x1;
   TVector3 x23 = x3 - x2;
   x12.SetZ(0.);
   x13.SetZ(0.);
   x23.SetZ(0.);
   Double_t x12mag = x12.Mag();
   x12 = x12.Unit();
   Double_t x13mag = x13.Mag();
   x13 = x13.Unit();
   Double_t x23mag = x23.Mag();
   x23 = x23.Unit();

   Double_t sinHalfPhi23 = x12.Cross(x13).Z();
   Double_t cosHalfPhi23 = 0.5 * (x13mag/x12mag + (1. - x23mag/x12mag)*(x12mag + x23mag)/x13mag);
   Double_t halfPhi23 = TMath::ATan2(sinHalfPhi23, cosHalfPhi23);

   Double_t r  = -0.5 * x23mag / sinHalfPhi23;

   TVector3 xc = 0.5 * (x2 + x3) + r * cosHalfPhi23 * x23.Cross(ez);

   if (dir == kIterBackward) r = -r;
   TMatrixD sv(5,1);
   sv(0,0) = 0;
   sv(1,0) = TMath::ATan2(r * (xc.Y() - x1.Y()), r * (xc.X() - x1.X()));
   sv(2,0) = GetPtoR() / r;
   sv(3,0) = 0;
   sv(4,0) = (x2.Z() - x3.Z()) / (r * 2 * halfPhi23);

   SetTo(sv, x1);
}

void THelicalTrack::CalcDapDa(Double_t  fid,
                              Double_t  dr,
                              Double_t  drp,
                              TMatrixD &F) const
{
   // ---------------------------------------------------
   // (2) Calculate @a'/@a = @a'/a = F_k-1
   // ---------------------------------------------------
   //        a' = (dr', fi0', cpa', dz', tnl')
   //        a  = (dr , fi0 , cpa , dz , tnl )
   //
                                                                                
   // @drho'/@a
   Double_t cpa   = fKappa;
   Double_t tnl   = fTanL;
   Double_t csfd  = TMath::Cos(fid);
   Double_t snfd  = TMath::Sin(fid);
   Double_t r     = fAlpha / cpa;
   Double_t rdr   = r + dr;
   Double_t rcpar = r / cpa;
   Double_t rdrpr = 1. / (r + drp);
                                                                                
   F(0,0) = csfd;
   F(0,1) = rdr*snfd;
   F(0,2) = rcpar*(1.-csfd);
   F(0,3) = 0.;
   F(0,4) = 0.;
                                                                                
   // @phi0'/@a
   F(1,0) = -rdrpr*snfd;
   F(1,1) =  rdr*rdrpr*csfd;
   F(1,2) =  rcpar*rdrpr*snfd;
   F(1,3) =  0.;
   F(1,4) =  0.;
                                                                                
   // @kappa'/@a
   F(2,0) = 0.;
   F(2,1) = 0.;
   F(2,2) = 1.;
   F(2,3) = 0.;
   F(2,4) = 0.;
                                                                                
   // @dz'/@a
   F(3,0) =  r*rdrpr*tnl*snfd;
   F(3,1) =  r*tnl*(1.-rdr*rdrpr*csfd);
   F(3,2) =  rcpar*tnl*(fid-r*rdrpr*snfd);
   F(3,3) =  1.;
   F(3,4) = -r*fid;
                                                                                
   // @tanl'/@a
   F(4,0) = 0.;
   F(4,1) = 0.;
   F(4,2) = 0.;
   F(4,3) = 0.;
   F(4,4) = 1.;
}

