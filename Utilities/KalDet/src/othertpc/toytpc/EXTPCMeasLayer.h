#ifndef EXTPCMEASLAYER_H
#define EXTPCMEASLAYER_H
//*************************************************************************
//* ===================
//*  EXTPCMeasLayer Class
//* ===================
//*
//* (Description)
//*   Sample measurement layer class used by EXTPCHit.
//* (Requires)
//* (Provides)
//*     class EXTPCMeasLayer
//* (Update Recored)
//*   2003/09/30  Y.Nakashima       Original version.
//*
//*************************************************************************
//
#include "TVector3.h"
#include "kaltest/TKalMatrix.h"
#include "kaltest/TPlane.h"
#include "kaltest/KalTrackDim.h"
#include "kaldet/EXVMeasLayer.h"

class TVTrackHit;

class EXTPCMeasLayer : public EXVMeasLayer, public TPlane {
public:
   // Ctors and Dtor

   EXTPCMeasLayer(TMaterial &min,
                  TMaterial &mout,
                  Double_t   y0,
                  Double_t   lhalf,
                  Double_t   sigmax0,
                  Double_t   sigmax1,
                  Double_t   sigmaz0,
                  Double_t   sigmaz1,
                  Bool_t     type = EXVMeasLayer::kActive);
   EXTPCMeasLayer(TMaterial &min,
                  TMaterial &mout,
                  Double_t   y0,
                  Double_t   xmin,
                  Double_t   ymax,
                  Double_t   lhalf,
                  Double_t   sigmax0,
                  Double_t   sigmax1,
                  Double_t   sigmaz0,
                  Double_t   sigmaz1,
                  Bool_t     type = EXVMeasLayer::kDummy,
                  Int_t      layer = -1,
                  Int_t      module = 0);
   virtual ~EXTPCMeasLayer();

   // Getters and Setters

   inline  Int_t GetModuleID() const { return fModule; }
   inline  Int_t GetLayerID () const { return fLayer;  }

   // Parrent's pure virtuals that must be implemented

   virtual TKalMatrix XvToMv    (const TVTrackHit &ht,
                                 const TVector3   &xv)   const;
   virtual TKalMatrix XvToMv    (const TVector3   &xv,
                                       Int_t       side) const;
   virtual TVector3   HitToXv   (const TVTrackHit &ht)   const;
   virtual void       CalcDhDa  (const TVTrackHit &ht,
                                 const TVector3   &xv,
                                 const TKalMatrix &dxphiada,
                                       TKalMatrix &H)    const;
   virtual void       ProcessHit(const TVector3   &xx,
                                       TObjArray  &hits);

   virtual Double_t GetSortingPolicy () const;

   Double_t GetSigmaX(Double_t z) const;
   Double_t GetSigmaZ(Double_t z) const;

private:
   Double_t fXMin;    // minimum phi
   Double_t fXMax;    // maximum phi
   Double_t fSigmaX0;   // xy resolution
   Double_t fSigmaX1;   // xy resolution
   Double_t fSigmaZ0;   // z  resolution
   Double_t fSigmaZ1;   // z  resolution
   Int_t    fModule;    // module number
   Int_t    fLayer;     // layer number

   ClassDef(EXTPCMeasLayer,1)   // Sample measurement layer class
};

#endif
