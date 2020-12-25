#ifndef ILDPLANARSTRIPHIT_H
#define ILDPLANARSTRIPHIT_H

/** ILDPlanarStripHit: User defined KalTest hit class using u coordinate, which provides coordinate vector as defined by the MeasLayer 
 *  
 * @author S.Aplin DESY
 */

#include "kaltest/KalTrackDim.h"

#include "ILDVTrackHit.h"


#define ILDPlanarStripHit_DIM 1 // set to 2 if one want to debug strip hits by using the 2nd dimention

class ILDPlanarStripHit : public ILDVTrackHit {
  
public:
  
  /** Constructor Taking a single coordinate and associated measurement layer, with bfield */
  ILDPlanarStripHit(const TVMeasLayer &ms,
               Double_t       *x,
               Double_t       *dx,
               Double_t        bfield,
               edm4hep::ConstTrackerHit trkhit) 
  : ILDVTrackHit(ms, x, dx, bfield, ILDPlanarStripHit_DIM,trkhit)
  { /* no op */ } 
  
  // TVTrackHit's pure virtuals that must be implemented
  
  /** Global to Local coordinates */
  virtual TKalMatrix XvToMv (const TVector3 &xv, Double_t t0) const;
  
  /** Print Debug information */
  virtual void       DebugPrint(Option_t *opt = "")           const;
  
  
private:
  
  
};
#endif
