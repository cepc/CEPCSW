#ifndef DETSEGMENTATION_GRIDDRIFTCHAMBER_H
#define DETSEGMENTATION_GRIDDRIFTCHAMBER_H

#include "DDSegmentation/Segmentation.h"

#include "TVector3.h"
#include <cmath>
#include <iostream>

/** GridDriftChamber Detector/DetSegmentation/DetSegmentation/GridDriftChamber.h GridDriftChamber.h
 *
 *  Segmentation for drift chamber.
 *
 *  @author    nalipour
 */

namespace dd4hep {
namespace DDSegmentation {
class GridDriftChamber : public Segmentation {
public:
  /// default constructor using an arbitrary type
  GridDriftChamber(const std::string& aCellEncoding);
  /// Default constructor used by derived classes passing an existing decoder
  GridDriftChamber(const BitFieldCoder* decoder);
  /// destructor
  virtual ~GridDriftChamber() = default;

  virtual Vector3D position(const CellID& aCellID) const;
  virtual CellID cellID(const Vector3D& aLocalPosition, const Vector3D& aGlobalPosition,
                        const VolumeID& aVolumeID) const;

//  inline double innerRadius() const { return m_innerRadius; }
//  inline double detectorLength() const { return m_detectorLength; }
  inline double offsetPhi() const { return m_offsetPhi; }
  inline double delta_phi() const{ return m_delta_phi; }
  inline const std::string& fieldNamePhi() const { return m_phiID; }
  // Setters

//  inline void setGeomParams(int layer, double sizePhi) {
//    layer_params[layer] = {sizePhi};
// }

//  void updateParams(int layer) const {
//    auto it_end = layer_params.cend();
//    --it_end;
//    double size = it_end->second[0];
//    double radius = it_end->second[1];
//    double eps = it_end->second[2];
//
//    auto map_it = layer_params.find(layer);
//    if (map_it != layer_params.cend()) {
//      size = map_it->second[0];
//      radius = map_it->second[1];
//      eps = map_it->second[2];
//    }
//
//    _currentGridSizePhi = size;
//    _currentRadius = radius;
//    m_epsilon = eps;
//  }

  inline double phiFromXY(const Vector3D& aposition) const {
    return std::atan2(aposition.Y, aposition.X) + M_PI ;
  }

//  inline int returnLayer(double x, double y) const {
//  // Hit R position
//    double R = std::sqrt(x * x + y * y);
//  // Layer
//    int layer = int((R - m_innerRadius) / m_cellSize);
//    return layer;
//  }

protected:
  /* *** nalipour *** */
  double phi(const CellID& cID) const;


  double m_offsetPhi;
  double m_delta_phi;
  std::string m_phiID;

  // Current parameters of the layer: sizePhi
//  mutable double _currentGridSizePhi;  // current size Phi
//  mutable double _currentRadius;       // current size radius
//  mutable double m_epsilon;
};
}
}
#endif /* DETSEGMENTATION_GRIDDRIFTCHAMBER_H */
