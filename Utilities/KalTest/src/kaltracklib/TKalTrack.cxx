//*************************************************************************
//* =================
//*  TKalTrack Class
//* =================
//*
//* (Description)
//*   Track class for Kalman filter
//* (Requires)
//*     TVKalSystem
//* (Provides)
//*     class TKalTrack
//* (Update Recored)
//*   2003/09/30  Y.Nakashima   Original version.
//*   2005/02/23  A.Yamaguchi   Added a new data member, fMass.
//*   2005/08/25  K.Fujii       Added drawable attribute.
//*   2005/08/26  K.Fujii       Removed drawable attribute.
//*
//*************************************************************************
                                                                                
#include "TKalTrackState.h"    // from KalTrackLib
#include "TKalTrackSite.h"     // from KalTrackLib
#include "TKalTrack.h"         // from KalTrackLib
#include <iostream>            // from STL

using namespace std;
#if __GNUC__ < 4 && !defined(__STRICT_ANSI__)
#else
const Double_t TKalTrack::kMpi = 0.13957018; // pion mass [GeV]
#endif


//_________________________________________________________________________
//  ------------------------------
//   TKalTrack: Kalman rack class
//  ------------------------------

ClassImp(TKalTrack)
                                                                                
//_________________________________________________________________________
//  ----------------------------------
//   Ctor
//  ----------------------------------
TKalTrack::TKalTrack(Int_t n)
          :TVKalSystem(n), fMass(kMpi)
{
}

//_________________________________________________________________________
//  ----------------------------------
//   Utility Method
//  ----------------------------------
//_________________________________________________________________________
// -----------------
//  FitToHelix
// -----------------
//    chi^2-fits hits belonging to this track to a single helix.
//
Double_t TKalTrack::FitToHelix(TKalTrackState &a, TKalMatrix &C, Int_t &ndf)
{
   // Define static constants...
 
   static const Double_t kChi2Dum = 1.e20;
   static const Double_t kChi2Tol = 1.e-8;
   static const Double_t kLambda  = 1.;
   static const Double_t kLincr   = 10.;
   static const Double_t kLdecr   = 0.1;
   static const Int_t    kLoopMax = 100;

   // Initialize return values...

   TKalTrackState abest(a);
   Double_t   chi2best = kChi2Dum;
   Double_t   lambda   = kLambda;
   Int_t      nloops   = 0;

   ndf = 0;
   Int_t      nsites   = 0;
   Double_t   chi2     = 0.;

   // Prepare som matrices

   TIter      next(this);
   TKalTrackSite &site0 = *(TKalTrackSite *)next();

   Int_t mdim = site0.GetDimension();
   Int_t sdim = site0.GetCurState().GetDimension();

   TKalMatrix    dchi2dabest(1   , sdim);
   TKalMatrix    dchi2da    (1   , sdim);
   TKalMatrix    d2chi2dada (sdim, sdim);
   TKalMatrix    d2chi2best (sdim, sdim);
   TKalMatrix    curh       (mdim, 1   );
   TKalMatrix    curH       (mdim, sdim);
   TKalMatrix    curHt      (sdim, mdim);
   TKalMatrix    curVinv    (mdim, mdim);
   TKalMatrix    curResVec  (mdim, 1   );

   // Minimization loop starts here

   while (1) {
       if (nloops > kLoopMax) {
          cerr << "TKalTrack::FitToHelix >>>>>>>>>>>>>>"
               << " Loop count limit reached. nloops = " << nloops << endl;
          a = abest;
          chi2  = chi2best;
          break;
       }
       nloops++;

       d2chi2dada.Zero();
       dchi2da.Zero();
       chi2 = 0;

       // Loop over hits and accumulate chi2

       next.Reset();
       TKalTrackSite *sitePtr = 0;
       nsites  = 0;
       while ((sitePtr = (TKalTrackSite *)next())) {
           TKalTrackSite   &site = *sitePtr;

           if (site.IsLocked())                      continue;
           if (!site.CalcExpectedMeasVec  (a, curh)) continue;
           if (!site.CalcMeasVecDerivative(a, curH)) continue;
           nsites++;	// site accepted

           //dchi2/da
           curVinv   = TKalMatrix(TKalMatrix::kInverted,site.GetMeasNoiseMat());
           curResVec = site.GetMeasVec() - curh;

           TKalMatrix curResVecT(TKalMatrix::kTransposed, curResVec);
           dchi2da += (curResVecT * curVinv * curH);

           // accumulate chi2
           Double_t delchi2 = (curResVecT * curVinv * curResVec)(0,0);
           chi2 += delchi2;

           // inverse covariance matrix
           TKalMatrix curHt(TKalMatrix::kTransposed, curH);
           d2chi2dada += (curHt * curVinv * curH);
       }

       //
       // if (chi2best - chi2) < kChi2Tol, break while loop.
       //

       if (TMath::Abs(chi2best - chi2) < kChi2Tol) {
           d2chi2best = d2chi2dada;
           break;
       }

       if (chi2 < chi2best) {
           // chi2 decreased. Save this step as the current best
           chi2best    = chi2;
           abest       = (TKalTrackState &)a;
           dchi2dabest = dchi2da;
           d2chi2best  = d2chi2dada;
           lambda     *= kLdecr;
       } else {
           // chi2 increased. Restore the current best
           dchi2da     = dchi2dabest;
           d2chi2dada  = d2chi2best;
           lambda     *= kLincr;
       }

       // Modify the 2nd derivative and try again

       for (Int_t i=0; i<sdim; i++) {
           d2chi2dada(i, i) *= (1 + lambda);
       }

       TKalMatrix d2chi2dadainv(TKalMatrix::kInverted, d2chi2dada);
       TKalMatrix dchi2daT     (TKalMatrix::kTransposed, dchi2da);
       a += (d2chi2dadainv * dchi2daT);
   }

   ndf  = nsites * mdim - sdim;
   C    = TKalMatrix(TKalMatrix::kInverted, d2chi2best);

   return chi2;
}
