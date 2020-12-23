#ifndef ILDPLANARHIT_H
#define ILDPLANARHIT_H

/** ILDPlanarHit: User defined KalTest hit class using u and v coordinates, which provides coordinate vector as defined by the MeasLayer 
 *
 * @author S.Aplin DESY
 */

#include "kaltest/KalTrackDim.h"

#include "ILDVTrackHit.h"

#define ILDPlanarHit_DIM 2

class ILDPlanarHit : public ILDVTrackHit {
  
public:
  
  /** Constructor Taking u and v coordinates and associated measurement layer, with bfield */
  ILDPlanarHit(const TVMeasLayer  &ms,
               Double_t           *x,
               Double_t           *dx,
               Double_t           bfield,
               edm4hep::ConstTrackerHit trkhit) 
  : ILDVTrackHit(ms, x, dx, bfield, ILDPlanarHit_DIM,trkhit)
  { /* no op */ } 
  
  // TVTrackHit's pure virtuals that must be implemented
  
  /** Global to Local coordinates */
  virtual TKalMatrix XvToMv (const TVector3 &xv, Double_t t0) const;
  
  /** Print Debug information */
  virtual void       DebugPrint(Option_t *opt = "")           const;
  
  
private:
  
  
};
#endif
