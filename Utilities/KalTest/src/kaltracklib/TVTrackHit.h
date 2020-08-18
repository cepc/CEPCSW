#ifndef TVTRACKHIT_H
#define TVTRACKHIT_H

//*************************************************************************
//* ====================
//*  TVTrackHit Class
//* ====================
//*
//* (Description)
//*   Abstract base class to store single hit information.
//* (Requires)
//* (Provides)
//*     class TVTrackHit
//* (Update Recored)
//*   2003/09/30  Y.Nakashima       Original version.
//*   2005/08/11  K.Fujii           Removed fXX and its getter and setter.
//*
//*************************************************************************

#include "TVector3.h"      // from ROOT
#include "TKalMatrix.h"    // from KalLib
#include "KalTrackDim.h"   // from KalTrackLib
#include "TVMeasLayer.h"   // from KalTrackLib

class TVTrackHit : public TKalMatrix {
public:
   TVTrackHit(Int_t m = kMdim);

   TVTrackHit(const TVMeasLayer &ms, Double_t *x, Double_t *dx, 
              Double_t b = 30., Int_t m = kMdim);
   TVTrackHit(const TVTrackHit &hit);

   virtual ~TVTrackHit();

   inline virtual Double_t GetX (Int_t i) const { return (*this)(i,0);      }
   inline virtual Double_t GetDX(Int_t i) const { return (*this)(i,1);      }
   inline virtual Int_t    GetDimension() const { return fDim;              }
   inline virtual Double_t GetBfield()    const { return fBfield;           }

   inline virtual const TVMeasLayer & GetMeasLayer() const 
                                                { return *fMeasLayerPtr;    }

   virtual TKalMatrix XvToMv  (const TVector3 &xv, Double_t t0)  const = 0;

  //   virtual void       DebugPrint(Option_t *opt = "")        const = 0;

private:
   Int_t         fDim;            // dimension of coordinate space
   Double_t      fBfield;         // B field
   TVMeasLayer  *fMeasLayerPtr;   // pointer to measurement layer

   ClassDef(TVTrackHit,1)      // Sample hit class
};

#endif
