#ifndef __EXEVENTGEN__
#define __EXEVENTGEN__

#include "kaltest/TKalDetCradle.h"
#include "kaltest/THelicalTrack.h"
#include "TMath.h"

class EXEventGen {
public:
   EXEventGen(TKalDetCradle const &cradle, TObjArray &kalhits)
             : fCradlePtr(&cradle), fHitBufPtr(&kalhits) {}
   virtual ~EXEventGen() {}

   THelicalTrack GenerateHelix(Double_t pt,
                               Double_t cosmin,
                               Double_t cosmax,
                               Double_t phimin=0.,
                               Double_t phimax=2*TMath::Pi(),
                               TVector3 xv0=TVector3(0.,0.,0.));
   void          Swim(THelicalTrack &heltrk);

   static void     SetT0(Double_t t0) { fgT0 = t0;   }
   static Double_t GetT0()            { return fgT0; }

private:
   TKalDetCradle const *fCradlePtr;     // pointer to detector system
   TObjArray     *fHitBufPtr;     // pointer to hit array

   static Double_t  fgT0;         // t0

};

#endif
