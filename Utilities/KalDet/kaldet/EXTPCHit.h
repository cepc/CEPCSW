#ifndef LCTPC_EXTPCHIT_H
#define LCTPC_EXTPCHIT_H

#include "GearTPCCylinderHit.h"
#include <kaltest/TVMeasLayer.h>

/**
 * A backward compatibility class for GearTPCCylinderHit.
 * Do not use this in new code, but use GearTPCCylinderHit directly. 
 * This class extends the GearTPCCylinderHit by a side, which is never used anywhere.
 *
 * \deprecated EXTPCHit
 */

class EXTPCHit : public kaldet::GearTPCCylinderHit
{
  public:
  /// The default constructor. 
  EXTPCHit(Int_t m = kMdim);

  /// Constructor initialising the original hit as 3D coordinates.
  EXTPCHit(const TVMeasLayer &ms,
                 Double_t       *x,
                 Double_t       *dx,
                 Int_t           side,
                 Double_t        v,
           const TVector3       &xx,
                 Double_t        b,
	         Int_t           m = kMdim);

  /// Constructor initialising the original hit with a reference pointer.
  EXTPCHit(const TVMeasLayer &ms,
                 Double_t       *x,
                 Double_t       *dx,
                 Int_t           side,
                 Double_t        v,
           const void           *hitp,
                 Double_t        b,
                 Int_t           m = kMdim);

  /// The destructor.
  virtual ~EXTPCHit();

  /// Get the side value which has been set in the constructor.
  inline       Int_t     GetSide  () const { return fSide;   }

 private: 
  Int_t           fSide;    /// (-1, +1) = (-z side, +z side)

  //  ClassDef(EXTPCHit, 1)  // EXTPC hit class

};

#endif // LCTPC_EXTPCHIT_H
