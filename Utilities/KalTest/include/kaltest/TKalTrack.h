#ifndef TKALTRACK_H
#define TKALTRACK_H
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
//*   2003/09/30  Y.Nakashima	Original version.
//*   2005/02/23  A.Yamaguchi	Added a new data member fMass and its
//*                             getter and setter.
//*   2005/08/15  K.Fujii       Removed fDir and its getter and setter.
//*   2005/08/25  K.Fujii       Added Drawable attribute.
//*   2005/08/26  K.Fujii       Removed Drawable attribute.
//*
//*************************************************************************
                                                                                
#include "TVKalSystem.h"       // from KalLib
#include "TKalTrackState.h"    // from KalTrackLib

//_________________________________________________________________________
//  ------------------------------
//   TKalTrack: Kalman Track class
//  ------------------------------
                                                                                
class TKalTrack : public TVKalSystem {
public:
   TKalTrack(Int_t n = 1);
   ~TKalTrack() {}

   inline virtual void      SetMass(Double_t m)           { fMass = m;    }
   inline virtual Double_t  GetMass()             const   { return fMass; }

   Double_t FitToHelix(TKalTrackState &a, TKalMatrix &C, Int_t &ndf);

private:
   Double_t     fMass;        // mass [GeV]

#if __GNUC__ < 4 && !defined(__STRICT_ANSI__)
   static const Double_t kMpi = 0.13957018; //! pion mass [GeV]
#else
   static const Double_t kMpi;              //! pion mass [GeV]
#endif

   ClassDef(TKalTrack,1)  // Base class for Kalman Filter
};

#endif
