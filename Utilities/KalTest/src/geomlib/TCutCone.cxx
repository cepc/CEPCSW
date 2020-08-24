//*************************************************************************
//* ====================
//*  TCutCone Class
//* ====================
//*
//* (Description)
//*   A class to implement a hyperboloidal surface object.
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
#include <iostream>
#include "TCircle.h"
#include "TCutCone.h"
#include "TVTrack.h"

using namespace std;

#if __GNUC__ < 4 && !defined(__STRICT_ANSI__)
#else
const Double_t TCutCone::kTol = 1.e-5; // tolerance
#endif

//_____________________________________________________________________
//  -----------------------------------
//  TCutCone Class
//  -----------------------------------

ClassImp(TCutCone)

//_____________________________________________________________________
//  -----------------------------------
//  Calculate S
//  -----------------------------------
//
Double_t TCutCone::CalcS(const TVector3 &xx) const
{
   TVector3 xxc = xx - fXc;
   Double_t r   = xxc.Perp();
   Double_t s   = (r - xxc.Z()*fTanA) * (r + xxc.Z()* fTanA);
   return s;
}

//_____________________________________________________________________
//  -----------------------------------
//  Calculate (@S/@x)
//  -----------------------------------
//
TMatrixD TCutCone::CalcDSDx(const TVector3 &xx) const
{
   TVector3 xxc = xx - fXc;
   TMatrixD dsdx(1,3);
   dsdx(0,0) =  2.* xxc.X();
   dsdx(0,1) =  2.* xxc.Y();
   dsdx(0,2) = -2.* xxc.Z() * fTanA * fTanA;
   return dsdx;
}
