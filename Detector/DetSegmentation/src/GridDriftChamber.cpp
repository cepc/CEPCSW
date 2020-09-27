#include "DetSegmentation/GridDriftChamber.h"

namespace dd4hep {
namespace DDSegmentation {

/// default constructor using an encoding string
GridDriftChamber::GridDriftChamber(const std::string& cellEncoding) : Segmentation(cellEncoding) {
  // define type and description
  _type = "GridDriftChamber";
  _description = "Drift chamber segmentation in the global coordinates";

  registerIdentifier("identifier_phi", "Cell ID identifier for phi", m_phiID, "phi");
  registerParameter("delta_phi", "delta phi", m_delta_phi, 0., SegmentationParameter::LengthUnit);
}

GridDriftChamber::GridDriftChamber(const BitFieldCoder* decoder) : Segmentation(decoder) {
  // define type and description
  _type = "GridDriftChamber";
  _description = "Drift chamber segmentation in the global coordinates";

  registerIdentifier("identifier_phi", "Cell ID identifier for phi", m_phiID, "phi");
  registerParameter("delta_phi", "delta phi", m_delta_phi, 0., SegmentationParameter::LengthUnit);
}

Vector3D GridDriftChamber::position(const CellID& /*cID*/) const {
  Vector3D cellPosition = {0, 0, 0};
  return cellPosition;
}

CellID GridDriftChamber::cellID(const Vector3D& /*localPosition*/, const Vector3D& globalPosition,
                                const VolumeID& vID) const {

  CellID cID = vID;

  double phi_hit = phiFromXY(globalPosition);
  double posx = globalPosition.X;
  double posy = globalPosition.Y;

  int lphi = (int) (phi_hit/m_delta_phi);
  _decoder->set(cID, m_phiID, lphi);

//  std::cout << " myliu: "
//            << " x: " << posx
//            << " y: " << posy
////            << " pre: " << phi_pre
//            << " phi_hit: " << phi_hit
//            << " lphi: " << lphi
//            << std::endl;
  return cID;
}


REGISTER_SEGMENTATION(GridDriftChamber)
}
}
