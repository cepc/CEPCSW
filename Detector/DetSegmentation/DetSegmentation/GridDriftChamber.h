#ifndef DETSEGMENTATION_GRIDDRIFTCHAMBER_H
#define DETSEGMENTATION_GRIDDRIFTCHAMBER_H

#include "DDSegmentation/Segmentation.h"

#include "TVector3.h"
#include <cmath>
#include <iostream>
#include <map>

/** GridDriftChamber Detector/DetSegmentation/DetSegmentation/GridDriftChamber.h GridDriftChamber.h
 *
 *  Segmentation for drift chamber.
 *
 *  @author    nalipour
 */


typedef struct Layer
 {
   double layerphi;
   double R;
   double eps;
   double offset;
   Layer(){};
   Layer(double x, double y, double z, double k):layerphi(x),R(y),eps(z),offset(k){};
   bool operator < (const Layer &a) const
   {
      return layerphi < a.layerphi;
   }
 } LAYER;

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
  virtual double distanceTrackWire(const CellID& cID, const TVector3& hit_start/*, const TVector3& hit_end*/) const;

  double phi(const CellID& cID) const;
  inline double cell_Size() const { return m_cellSize; }
  inline double epsilon0() const { return m_epsilon0; }
  inline double detectorLength() const { return m_detectorLength; }
  inline const std::string& fieldNamePhi() const { return m_phiID; }
  // Setters

  inline double phiFromXY(const Vector3D& aposition) const {
    double hit_phi =  std::atan2(aposition.Y, aposition.X) ;
    if( hit_phi < 0 ) { hit_phi += 2 * M_PI; }
    return hit_phi;
  }

  inline void setGeomParams(int layer, double layerphi, double R, double eps, double offset) {
    layer_params.insert(std::pair<int,LAYER>(layer,LAYER(layerphi,R,eps,offset)));
   }

  inline void setWiresInLayer(int layer, int numWires)
  {
    double phi0;
    updateParams(layer);
    for (int i = 0; i<numWires; ++i) {

      double phi0 = m_offset;

      auto phi_start = _currentLayerphi * i + phi0;
      auto phi_end = phi_start + _currentLayerphi;

      TVector3 Wstart = returnWirePosition(phi_start, 1);
      TVector3 Wend = returnWirePosition(phi_end, -1);

      TVector3 Wmid = (Wstart+Wend)*(1/2.0);
      TVector3 Wdirection = (Wend - Wstart);

      m_wiresPositions[layer].push_back(std::make_pair(Wmid, Wdirection));
      }
  }

  inline auto returnAllWires() const { return m_wiresPositions; }

  inline TVector3 returnWirePosition(double angle, int sign) const {
    TVector3 w(0, 0, 0);
    w.SetX(_currentRadius * std::cos(angle));
    w.SetY(_currentRadius * std::sin(angle));
    w.SetZ(sign * m_detectorLength / 2.0);
    return w;
  }

  void updateParams(int layer)  const{
    auto it_end = layer_params.cend();
    --it_end;
    double layerphi = it_end->second.layerphi;
    double radius = it_end->second.R;
    double eps = it_end->second.eps;
    double offset = it_end->second.offset;

    auto map_it = layer_params.find(layer);
    if (map_it != layer_params.cend()) {
     layerphi = map_it->second.layerphi;
     radius = map_it->second.R;
     eps = map_it->second.eps;
     offset = map_it->second.offset;
    }
    _currentLayerphi = layerphi;
    _currentRadius = radius;
    m_epsilon = eps;
    m_offset = offset;
 }

 inline double returnAlpha() const {
   double alpha = 2 * std::asin(m_detectorLength * std::tan(m_epsilon0)/(2 * _currentRadius));
   return alpha;
 }

protected:

  std::map<int,LAYER> layer_params; // <layer, {layerphi, R, eps, offset}>
  std::map<int, std::vector<std::pair<TVector3, TVector3> >> m_wiresPositions; // < layer, vec<WireMidpoint, WireDirection> >

  double m_cellSize;
  double m_epsilon0;
  double m_detectorLength;
  std::string m_phiID;

  // Current parameters of the layer: sizePhi
  mutable double _currentLayerphi;
  mutable double _currentRadius;
  mutable double m_epsilon;
  mutable double m_offset;

};
}
}
#endif /* DETSEGMENTATION_GRIDDRIFTCHAMBER_H */
