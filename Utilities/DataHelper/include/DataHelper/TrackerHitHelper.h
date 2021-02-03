#ifndef TrackerHitHelper_H
#define TrackerHitHelper_H

#include "edm4hep/TrackerHit.h"
#include <array>

namespace CEPC{
  std::array<float, 6> GetCovMatrix(edm4hep::TrackerHit& hit);
  float                GetResolutionRPhi(edm4hep::TrackerHit& hit);
  float                GetResolutionZ(edm4hep::TrackerHit& hit);
}

#endif
