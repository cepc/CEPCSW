#ifndef GEARTPCCYLINDERHIT_H
#define GEARTPCCYLINDERHIT_H

#include <kaltest/KalTrackDim.h>
#include "GearTPCHit.h"
#include <kaltest/TVMeasLayer.h>

namespace kaldet{

/** The cylindrical implementation of the GearTPCHit.
 */
class GearTPCCylinderHit : public GearTPCHit {

public:
    /// KILLENB What does this constructor do? Best throw it out, it does not 
    /// properly initialise the class at all, does it?
  GearTPCCylinderHit(Int_t m = kMdim);

  /** Constructor to initialise the GearTPCHit using space point coordinates (TVector3) as original hit.
   */
  GearTPCCylinderHit(const TVMeasLayer &ms,
                 Double_t       *x,
                 Double_t       *dx,
           const TVector3       &xx,
                 Double_t        b,
                 Double_t        v,
                 Int_t           m = kMdim);

  /** Constructor using a pointer to the original hit as reference.
   */
  GearTPCCylinderHit(const TVMeasLayer &ms,
                 Double_t       *x,
                 Double_t       *dx,
           const void           *hitp,
                 Double_t        b,
                 Double_t        v,
                 Int_t           m = kMdim);

  /** The dectructor.
   */
  virtual ~GearTPCCylinderHit();

  /** Implementation of the space vector (xv) to measurement vector (mv) calculation
   *  for a cylindrical hit.
   */
  virtual TKalMatrix XvToMv(const TVector3 &xv, Double_t t0) const;
  
  /** Print some debug output to std err.
   */
  virtual void       DebugPrint(Option_t *opt = "")          const;
};

}//namespace kaldet

#endif //GEARTPCCYLINDERHIT_H
