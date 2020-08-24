//*************************************************************************
//* ======================
//*  TKalTrackState Class
//* ======================
//*
//* (Description)
//*   Track state vector class used in Kalman Filter.
//* (Requires)
//*     TVKalState
//* (Provides)
//*     class TKalTrackState
//* (Update Recored)
//*   2003/09/30  Y.Nakashima       Original version.
//*   2005/02/23  Y.Yamaguchi       Improved CalcProcessNoise().
//*   2005/08/14  K.Fujii           Removed CalcProcessNoise() and
//*                                 let TKalDetCradle::Transport() do its
//*                                 function.
//*   2010/04/06  K.Fujii           Modified MoveTo to allow a 1-dim hit,
//*                                 for which pivot is at the xpected hit.
//*
//*************************************************************************

#include "TKalDetCradle.h"      // from KalTrackLib
#include "TVKalDetector.h"      // from KalTrackLib
#include "TKalTrackState.h"     // from KalTrackLib
#include "TKalTrackSite.h"      // from KalTrackLib
#include "TKalTrack.h"          // from KalTrackLib

#include <iostream>             // from STL
#include <memory>               // from STL

using namespace std;

//_________________________________________________________________________
// -----------------------------------
//  Kalman track state vector
// -----------------------------------
//_________________________________________________________________________
// -----------------
//  Ctors and Dtor
// -----------------
//
TKalTrackState::TKalTrackState(Int_t p) 
           : TVKalState(p), fX0()
{
}

TKalTrackState::TKalTrackState(const TKalMatrix &sv, Int_t type, Int_t p) 
           : TVKalState(sv,type,p), fX0()
{
}

TKalTrackState::TKalTrackState(const TKalMatrix &sv, const TKalMatrix &c,
                                     Int_t type, Int_t p) 
           : TVKalState(sv,c,type,p), fX0()
{
}

TKalTrackState::TKalTrackState(const TKalMatrix &sv, const TVKalSite &site,
                                     Int_t type, Int_t p) 
           : TVKalState(sv,site,type,p), 
             fX0(((TKalTrackSite *)&site)->GetPivot())
{
}

TKalTrackState::TKalTrackState(const TKalMatrix &sv, const TKalMatrix &c,
                               const TVKalSite &site, Int_t type, Int_t p) 
           : TVKalState(sv,c,site,type,p),
             fX0(((TKalTrackSite *)&site)->GetPivot())
{
}

//_________________________________________________________________________
// ----------------------------------------------
//  Implementation of base-class pure virtuals
// ----------------------------------------------
//
TKalTrackState * TKalTrackState::MoveTo(TVKalSite  &to,
                                        TKalMatrix &F,
                                        TKalMatrix *QPtr) const
{	
   if (QPtr) {
      const TKalTrackSite &from   = static_cast<const TKalTrackSite &>(GetSite());
            TKalTrackSite &siteto = static_cast<TKalTrackSite &>(to);
            TKalDetCradle &det    = const_cast<TKalDetCradle &>
                                       (static_cast<const TKalDetCradle &>
                                          (from.GetHit().GetMeasLayer().GetParent()));
      Int_t sdim = GetDimension();
      TKalMatrix sv(sdim,1);
      det.Transport(from, siteto, sv, F, *QPtr); // siteto's pivot might be modified
      if (sdim == 6) {
         sv(5,0) = (*this)(5,0);
         F (5,5) = 1.;
      }
      return new TKalTrackState(sv, siteto, TVKalSite::kPredicted, sdim);
   } else {
      return 0;
   }
}

TKalTrackState & TKalTrackState::MoveTo(TVKalSite  &to,
                                        TKalMatrix &F,
                                        TKalMatrix &Q) const
{
   return *MoveTo(to, F, &Q);
}

void TKalTrackState::DebugPrint() const
{
   cerr << "          +-     -+   " << "+-" <<  endl
        << "          | drho  |   " << "| " << (*this)(0,0) << endl
        << "          | phi0  |   " << "| " << (*this)(1,0) << endl
        << " a      = | kappa | = " << "| " << (*this)(2,0) << endl
        << "          | dz    |   " << "| " << (*this)(3,0) << endl
        << "          | tanl  |   " << "| " << (*this)(4,0) << endl;
   if (GetDimension() == 6) {
      cerr 
        << "          | t0    |   " << "| " << (*this)(5,0) << endl;
   }
   cerr << "          +-     -+   " << "+-" << endl;
   cerr << "          +-" << endl 
        << " X0     = | " << fX0.X() << endl
        << "          | " << fX0.Y() << endl
        << "          | " << fX0.Z() << endl
        << "          +-" << endl;
   GetCovMat().DebugPrint(" covMat = ", 6);
}

//_________________________________________________________________________
// --------------------------------
//  Derived class methods
// --------------------------------
//
THelicalTrack TKalTrackState::GetHelix() const
{
   TKalMatrix a(5,1);
   for (Int_t i=0; i<5; i++) a(i,0) = (*this)(i,0);
   return THelicalTrack(a,fX0,((TKalTrackSite *)&GetSite())->GetBfield());
}

TStraightTrack TKalTrackState::GetLine() const
{
   TKalMatrix a(5,1);
   for (Int_t i=0; i<5; i++) a(i,0) = (*this)(i,0);
   return TStraightTrack(a,fX0,((TKalTrackSite *)&GetSite())->GetBfield());
}

TVTrack &TKalTrackState::CreateTrack() const
{
   TVTrack *tkp = 0;

   TKalMatrix a(5,1);
   for (Int_t i=0; i<5; i++) a(i,0) = (*this)(i,0);
   Double_t bfield = static_cast<const TKalTrackSite &>(GetSite()).GetBfield();

   if (bfield == 0.) tkp = new TStraightTrack(a,fX0);
   else              tkp = new THelicalTrack(a,fX0, bfield);

   return *tkp;
}

