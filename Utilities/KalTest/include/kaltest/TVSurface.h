#ifndef TVSURFACE_H
#define TVSURFACE_H
//*************************************************************************
//* ====================
//*  TVSurface Class
//* ====================
//*
//* (Description)
//*   This is the base class for various solids.
//* (Requires)
//*     TObject;
//* (Provides)
//*     class TVSurface
//* (Update Recored)
//*   2003/10/03  K.Fujii       Original version.
//*   2005/02/23  K.Fujii       Added new methods, Compare() and
//*                             GetSortingPolicy().
//*
//*   2011/06/17  D.Kamai       Added new method, GetOutwardNormal() 
//*                             
//*************************************************************************
//
#include "TObject.h"
#include "TMatrixD.h"
#include "TVector3.h"

class TVTrack;
//_____________________________________________________________________
//  -----------------------------------
//  Base Class for any surface
//  -----------------------------------

class TVSurface : public TObject {
public:

   virtual Int_t    CalcXingPointWith(const TVTrack  &hel,
                                            TVector3 &xx,
                                            Double_t &phi,
                                            Double_t  eps = 1.e-8) const;
   virtual Int_t    CalcXingPointWith(const TVTrack  &hel,
                                            TVector3 &xx,
                                            Double_t &phi,
                                            Int_t     mode,
   				            Double_t  eps = 1.e-8) const;

   virtual Double_t CalcS            (const TVector3 &xx) const = 0;
   virtual TMatrixD CalcDSDx         (const TVector3 &xx) const = 0;
   virtual Bool_t   IsOnSurface      (const TVector3 &xx) const = 0;
   virtual Bool_t   IsOutside        (const TVector3 &xx) const = 0;
   inline virtual TVector3 GetOutwardNormal (const TVector3 &xx) const;
  
   virtual Double_t GetSortingPolicy ()                   const = 0;

   virtual Int_t    Compare   (const TObject *obj) const;
   virtual Bool_t   IsSortable()                   const { return kTRUE; }
   
private:
 
   ClassDef(TVSurface,1)      // Base class for any surface
};

TVector3 TVSurface::GetOutwardNormal(const TVector3 &xx) const
{
  TMatrixD dsdx = CalcDSDx(xx);
  return TVector3(dsdx(0,0),dsdx(0,1),dsdx(0,2)).Unit();
}

#endif
