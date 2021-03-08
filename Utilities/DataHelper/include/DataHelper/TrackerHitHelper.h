#ifndef TrackerHitHelper_H
#define TrackerHitHelper_H

#include "edm4hep/TrackerHit.h"
#include <array>

namespace CEPC{
  std::array<float, 6> GetCovMatrix(edm4hep::TrackerHit& hit, bool useSpacePointerBuilderMethod = false);
  float                GetResolutionRPhi(edm4hep::TrackerHit& hit);
  float                GetResolutionZ(edm4hep::TrackerHit& hit);
  std::array<float, 6> ConvertToCovXYZ(float dU, float thetaU, float phiU, float dV, float thetaV, float phiV, bool useSpacePointBuilderMethod = false);
}

#endif
