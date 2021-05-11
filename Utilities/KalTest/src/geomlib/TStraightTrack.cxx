//*************************************************************************
//* ======================
//*  TStraightTrack Class
//* ======================
//*
//* (Description)
//*   A class to implement a straight track object.
//* (Requires)
//*     TVTrack, TCylinder, TCircle, TVector3, TMatrixD
//* (Provides)
//*     class TStraightTrack
//* (Update Recored)
//*   2003/10/24  K.Fujii       Original version.
//*
//*************************************************************************
//

#include <iostream>
#include "TStraightTrack.h"

using namespace std;

//_____________________________________________________________________
//  -----------------------------------
//  Straight Track Class
//  -----------------------------------
//
//_____________________________________________________________________
//  --------------
//  Ctors and Dtor
//  --------------
//

TStraightTrack::TStraightTrack(Double_t dr,
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

TStraightTrack::TStraightTrack(const TMatrixD &a,
                               const TVector3 &x0,
                                     Double_t  b)
             : TVTrack(a, x0, b)
{
}

//_____________________________________________________________________
//  ----------------
//  Utility methods
//  ----------------

void TStraightTrack::MoveTo(const TVector3 &xv0to,
                                  Double_t &t,
                                  TMatrixD *FPtr,
                                  TMatrixD *CPtr)
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

   Double_t csf0  = TMath::Cos(fi0);
   Double_t snf0  = TMath::Sqrt( TMath::Max(0.0, (1.0-csf0)*(1.0+csf0)) );
   if (fi0 > kPi) snf0 = -snf0;

            t     = -(xv - x0) * snf0 + (yv - y0) * csf0;

   Double_t drp   = dr - (xv - x0) * csf0 - (yv - y0) * snf0;
   Double_t dzp   = dz - (zv - z0) + t * tnl;

   TMatrixD av(5,1);
   av(0,0) = drp;
   av(1,0) = fi0;
   av(2,0) = cpa;
   av(3,0) = dzp;
   av(4,0) = tnl;

   TStraightTrack helto(av,xv0to);
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

   // @drho'/@a
   F(0,0) =  1;
   F(0,1) = -t;
   F(0,2) =  0;
   F(0,3) =  0;
   F(0,4) =  0;
                                                                                
   // @phi0'/@a
   F(1,0) =  0;
   F(1,1) =  1;
   F(1,2) =  0;
   F(1,3) =  0;
   F(1,4) =  0;
                                                                                
   // @kappa'/@a

   F(2,0) =  0;
   F(2,1) =  0;
   F(2,2) =  1;
   F(2,3) =  0;
   F(2,4) =  0;
   
   // @dz'/@a
   F(3,0) =  0;
   F(3,1) =  (drp - dr)* tnl;
   F(3,2) =  0;
   F(3,3) =  1;
   F(3,4) =  t;
    
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

TVector3 TStraightTrack::CalcXAt(Double_t t) const
{
   Double_t csf0 = TMath::Cos(fPhi0);
   Double_t snf0 = TMath::Sin(fPhi0);

   Double_t x    = fX0.X() + fDrho * csf0 - t * snf0;
   Double_t y    = fX0.Y() + fDrho * snf0 + t * csf0;
   Double_t z    = fX0.Z() + fDz          + t * fTanL;

   return TVector3(x,y,z);
}

TMatrixD TStraightTrack::CalcDxDa(Double_t t) const
{
   Double_t fi0   = fPhi0;

   Double_t snf0 = TMath::Sin(fi0);
   Double_t csf0 = TMath::Cos(fi0);

   TMatrixD dxda(3,5);
   // @x/@a
   dxda(0,0) =  csf0;
   dxda(0,1) = - fDrho * snf0 - t * csf0;
   dxda(0,2) = 0;
   dxda(0,3) = 0;
   dxda(0,4) = 0;

   // @y/@a
   dxda(1,0) = snf0;
   dxda(1,1) = fDrho * csf0 - t * snf0;
   dxda(1,2) = 0;
   dxda(1,3) = 0;
   dxda(1,4) = 0;

   // @z/@a
   dxda(2,0) = 0;
   dxda(2,1) = 0;
   dxda(2,2) = 0;
   dxda(2,3) = 1;
   dxda(2,4) = t;

   return dxda;
}

TMatrixD TStraightTrack::CalcDxDphi(Double_t t) const
{
   Double_t snf0 = TMath::Sin(fPhi0);
   Double_t csf0 = TMath::Cos(fPhi0);

   TMatrixD dxdphi(3,1);
   dxdphi(0,0) = -snf0;
   dxdphi(1,0) =  csf0;
   dxdphi(2,0) =  fTanL;

   return dxdphi;
}

void TStraightTrack::CalcDapDa(Double_t fid,
                               Double_t dr,
                               Double_t drp,
                               TMatrixD &F) const
{
   // @drho'/@a
   F(0,0) =  1;
   F(0,1) = -fid;
   F(0,2) =  0;
   F(0,3) =  0;
   F(0,4) =  0;
                                                                                
   // @phi0'/@a
   F(1,0) =  0;
   F(1,1) =  1;
   F(1,2) =  0;
   F(1,3) =  0;
   F(1,4) =  0;
                                                                                
   // @kappa'/@a
                                                                                
   F(2,0) =  0;
   F(2,1) =  0;
   F(2,2) =  1;
   F(2,3) =  0;
   F(2,4) =  0;
                                                                                
   // @dz'/@a
   F(3,0) =  0;
   F(3,1) =  (drp - dr)* fTanL;
   F(3,2) =  0;
   F(3,3) =  1;
   F(3,4) =  fid;
                                                                                
   // @tanl'/@a
   F(4,0) = 0;
   F(4,1) = 0;
   F(4,2) = 0;
   F(4,3) = 0;
   F(4,4) = 1;
} 

