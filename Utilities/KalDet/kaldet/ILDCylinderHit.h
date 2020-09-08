#ifndef ILDCYLINDERHIT_H
#define ILDCYLINDERHIT_H

/** ILDCylinderHit: User defined KalTest hit class using R and Rphi coordinates, which provides coordinate vector as defined by the MeasLayer 
 *
 * @author S.Aplin DESY
 */

#include "kaltest/KalTrackDim.h"
#include "ILDVTrackHit.h"


class ILDCylinderHit : public ILDVTrackHit {
  
public:
  
  
  /** Constructor Taking R and Rphi coordinates and associated measurement layer, with bfield */
  ILDCylinderHit(const TVMeasLayer &ms, Double_t *x, Double_t *dx, 
                 Double_t bfield, edm4hep::ConstTrackerHit trkhit ) 
  : ILDVTrackHit(ms, x, dx, bfield, 2, trkhit)
  { /* no op */ } 
    
  
  // TVTrackHit's pure virtuals that must be implemented
  
  /** Global to Local coordinates */
  virtual TKalMatrix XvToMv(const TVector3 &xv, Double_t t0) const;
  
  /** Print Debug information */
  virtual void       DebugPrint(Option_t *opt = "")         const;
  
  
private:
  
  
};
#endif
