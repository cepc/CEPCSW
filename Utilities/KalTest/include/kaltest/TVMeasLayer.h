#ifndef TVMEASLAYER_H
#define TVMEASLAYER_H
//*************************************************************************
//* ====================
//*  TVMeasLayer Class
//* ====================
//*
//* (Description)
//*   Measurement layer interface class.
//* (Requires)
//* (Provides)
//*     class TVMeasLayer
//* (Update Recored)
//*   2003/09/30  Y.Nakashima       Original version.
//*   2005/02/23  A.Yamaguchi       Added new data members, fFwdX0Inv,
//*                                 fBwdX0Inv and fIndex, and their
//*                                 corresponding getters and setters.
//*                                 Added a new method, GetX0Inv().
//*   2005/0X/XX  A.Yamaguchi       Replaced fFwdX0Inv, fBwdX0Inv, and
//*                                 their getters and setters by
//*                                 fMaterialOutPtr, fMaterialInPtr, and
//*                                 their getters and setters.
//*   2005/08/15  K.Fujii           Added fIsActive and IsActive().
//*   2011/12/03  S.Aplin           Added new member: name 
//*                                 default value set to "TVMeasLayer"
//*                                 and corresponding member function
//*                                 TString GetName()  
//*
//*************************************************************************

#include "TVector3.h"       // from ROOT
#include "TMaterial.h"      // from ROOT
#include "TAttElement.h"    // from Utils
#include "TKalMatrix.h"     // from KalLib
#include "KalTrackDim.h"    // from KalTrackLib

class TVTrack;
class TVTrackHit;

class TVMeasLayer : public TAttElement {
public:
   // Ctors and Dtor

   TVMeasLayer(TMaterial &matIn, 
               TMaterial &matOut,
               Bool_t     isactive = kTRUE,
               const Char_t    *name = "TVMeasLayer");
   virtual ~TVMeasLayer() {}

   // Utiliy Methods

   virtual TKalMatrix XvToMv   (const TVTrackHit &ht,
                                const TVector3   &xv) const = 0;
   virtual TVector3   HitToXv  (const TVTrackHit &ht) const = 0;
   virtual void       CalcDhDa (const TVTrackHit &ht,
                                const TVector3   &xv,
                                const TKalMatrix &dxphiada,
                                      TKalMatrix &H)  const = 0;

   inline virtual TMaterial &GetMaterial(Bool_t isoutgoing) const
              { return isoutgoing ? *fMaterialOutPtr : *fMaterialInPtr; }
     
   inline  Int_t      GetIndex() const  { return fIndex;    }
   inline  void       SetIndex(Int_t i) { fIndex = i;       }    
   inline  Bool_t     IsActive() const  { return fIsActive; }

   virtual Double_t   GetEnergyLoss (      Bool_t    isoutgoing,
                                     const TVTrack  &hel,
                                           Double_t  df) const;
   virtual void       CalcQms       (      Bool_t    isoutgoing,
                                     const TVTrack  &hel,
                                           Double_t  df,
                                           TKalMatrix &Qms) const;

  inline TString       GetName() const { return fname;    }
  
private:
   TMaterial     *fMaterialInPtr;   // pointer of inner Material
   TMaterial     *fMaterialOutPtr;  // pointer of outer Material
   Int_t          fIndex;           // index in TKalDetCradle
   Bool_t         fIsActive;        // flag to tell layer is active or not
  const Char_t   *fname;
   ClassDef(TVMeasLayer,1)      // Measurement layer interface class
};

#endif
