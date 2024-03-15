#ifndef TrackerHitHelper_H
#define TrackerHitHelper_H

#include "edm4hep/TrackerHit.h"
#include "edm4hep/SimTrackerHit.h"
#include "edm4hep/MCRecoTrackerAssociationCollection.h"
#include "DDSegmentation/Segmentation.h"
#include "DetSegmentation/GridDriftChamber.h"
#include <array>

//namespace dd4hep {
//         class Detector;
//         namespace DDSegmentation{
//             class GridDriftChamber;
//         }
//}

namespace CEPC{
  std::array<float, 6> GetCovMatrix(edm4hep::TrackerHit& hit, bool useSpacePointerBuilderMethod = false);
  float                GetResolutionRPhi(edm4hep::TrackerHit& hit);
  float                GetResolutionZ(edm4hep::TrackerHit& hit);
  std::array<float, 6> ConvertToCovXYZ(float dU, float thetaU, float phiU, float dV, float thetaV, float phiV, bool useSpacePointBuilderMethod = false);
  const edm4hep::SimTrackerHit getAssoClosestSimTrackerHit(
          const edm4hep::MCRecoTrackerAssociationCollection* assoHits,
          const edm4hep::TrackerHit trackerHit,
          const dd4hep::DDSegmentation::GridDriftChamber* segmentation,
          int docaMehtod);
}

#endif
