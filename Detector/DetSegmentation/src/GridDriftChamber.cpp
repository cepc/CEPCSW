#include "DetSegmentation/GridDriftChamber.h"
#include <map>

namespace dd4hep {
namespace DDSegmentation {

/// default constructor using an encoding string
GridDriftChamber::GridDriftChamber(const std::string& cellEncoding) : Segmentation(cellEncoding) {
  // define type and description
  _type = "GridDriftChamber";
  _description = "Drift chamber segmentation in the global coordinates";

  registerParameter("cell_size", "cell size", m_cellSize, 0., SegmentationParameter::LengthUnit);
  registerParameter("detector_length", "Length of the wire", m_detectorLength, 1., SegmentationParameter::LengthUnit);
  registerIdentifier("identifier_phi", "Cell ID identifier for phi", m_phiID, "cellID");
}

GridDriftChamber::GridDriftChamber(const BitFieldCoder* decoder) : Segmentation(decoder) {

  _type = "GridDriftChamber";
  _description = "Drift chamber segmentation in the global coordinates";

  registerParameter("cell_size", "cell size", m_cellSize, 1., SegmentationParameter::LengthUnit);
  registerParameter("epsilon0", "epsilon", m_epsilon0, 0., SegmentationParameter::AngleUnit, true);
  registerParameter("detector_length", "Length of the wire", m_detectorLength, 1., SegmentationParameter::LengthUnit);
  registerIdentifier("identifier_phi", "Cell ID identifier for phi", m_phiID, "cellID");
}

Vector3D GridDriftChamber::position(const CellID& /*cID*/) const {
  Vector3D cellPosition = {0, 0, 0};
  return cellPosition;
}


CellID GridDriftChamber::cellID(const Vector3D& /*localPosition*/, const Vector3D& globalPosition,
                                const VolumeID& vID) const {

  CellID cID = vID;
  unsigned int layerID = _decoder->get(vID, "layer");
  updateParams(layerID);

  double phi_hit = phiFromXY(globalPosition);
  double posx = globalPosition.X;
  double posy = globalPosition.Y;
  double offsetphi= m_offset;
  int _lphi;

  if(phi_hit >= offsetphi) {
    _lphi = (int) ((phi_hit - offsetphi)/ _currentLayerphi);
  }
  else {
    _lphi = (int) ((phi_hit - offsetphi + 2 * M_PI)/ _currentLayerphi);
  }

  int lphi = _lphi;
  _decoder->set(cID, m_phiID, lphi);


std::cout << "#######################################: " 
          <<  " offset : " << m_offset
          << " offsetphi: " << offsetphi
          << " layerID: " << layerID
          << " r: " << _currentRadius
          << " layerphi: " << _currentLayerphi
          << std::endl;

  return cID;
}

double GridDriftChamber::phi(const CellID& cID) const {
  CellID phiValue = _decoder->get(cID, m_phiID);
  return binToPosition(phiValue, _currentLayerphi, m_offset);
}

void GridDriftChamber::cellposition(const CellID& cID, TVector3& Wstart,
                                    TVector3& Wend) const {

  auto layerIndex = _decoder->get(cID, "layer");
  updateParams(layerIndex);

  double phi_start = phi(cID);
  double phi_end = phi_start + returnAlpha();

  Wstart = returnWirePosition(phi_start, -1);
  Wend = returnWirePosition(phi_end, 1);
}



double GridDriftChamber::distanceTrackWire(const CellID& cID, const TVector3& hit_start,
                                           const TVector3& hit_end) const {

//  auto layerIndex = _decoder->get(cID, "layer");
//  updateParams(layerIndex);
//
//  double phi_start = phi(cID);
//  double phi_end = phi_start + returnAlpha();

//  TVector3 Wstart = returnWirePosition(phi_start, -1); // The default centimeter unit in DD4hep
//  TVector3 Wend = returnWirePosition(phi_end, 1);   // The default centimeter unit in DD4hep
  TVector3 Wstart = {0,0,0};
  TVector3 Wend = {0,0,0};
  cellposition(cID,Wstart,Wend);

  TVector3 a = hit_end - hit_start;
  TVector3 b = Wend - Wstart;
  TVector3 c = Wstart - hit_start;

  double num = std::abs(c.Dot(a.Cross(b)));
  double denum = (a.Cross(b)).Mag();
//  double num = (b.Cross(c)).Mag();
//  double denum = b.Mag();

  double DCA = 0;

   if (denum) {
    DCA = num / denum;
  }

  return DCA;
}


REGISTER_SEGMENTATION(GridDriftChamber)
}
}
