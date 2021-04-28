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
   Layer(double x, double y, double z, double k):layerphi(x),R(y),eps(z),offset(k){}
   bool operator < (const Layer &a) const
   {
      return layerphi < a.layerphi;
   }
 } LAYER;

typedef struct CID
 {
   int chamberID;
   int layerID;
//   CID(){}
   CID(int i, int j): chamberID(i),layerID(j){}
   // the operator < defines the operation used in map
   friend bool operator < (const CID &c1, const CID &c2);
 } vID;

inline bool operator < (const struct CID &c1, const struct CID &c2) {
    return c1.chamberID < c2.chamberID || (c1.chamberID == c2.chamberID && c1.layerID < c2.layerID);
}

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
//  virtual int LayerID(const Vector3D& aGlobalPosition) const;
  virtual CellID cellID(const Vector3D& aLocalPosition, const Vector3D& aGlobalPosition,
                        const VolumeID& aVolumeID) const;
  virtual double distanceTrackWire(const CellID& cID, const TVector3& hit_start, const TVector3& hit_end) const;
  virtual void cellposition(const CellID& cID, TVector3& Wstart, TVector3& Wend) const;

//  double phi(const CellID& cID) const;
  inline double cell_Size() const { return m_cellSize; }
  inline double epsilon0() const { return m_epsilon0; }
  inline double detectorLength() const { return m_detectorLength; }
  inline double DC_inner_rbegin() const { return m_DC_inner_rbegin; }
  inline double DC_inner_rend() const { return m_DC_inner_rend; }
  inline double DC_outer_rbegin() const { return m_DC_outer_rbegin; }
  inline double DC_outer_rend() const { return m_DC_outer_rend; }
  inline double DC_inner_rmin() const { return m_DC_inner_rmin; }
  inline double DC_inner_rmax() const { return m_DC_inner_rmax; }
  inline double DC_outer_rmin() const { return m_DC_outer_rmin; }
  inline double DC_outer_rmax() const { return m_DC_outer_rmax; }
  inline double safe_distance() const { return m_safe_distance; }
  inline double layer_width() const { return m_layer_width; }
  inline int DC_inner_layer_number() const { return m_DC_inner_layer_number; }
  inline int DC_outer_layer_number() const { return m_DC_outer_layer_number; }
  inline const std::string& fieldNamePhi() const { return m_phiID; }
  inline const std::string& Layerid() const { return layer_id; }
  // Setters

  inline double phiFromXY(const Vector3D& aposition) const {
    double hit_phi =  std::atan2(aposition.Y, aposition.X) ;
    if( hit_phi < 0 ) { hit_phi += 2 * M_PI; }
    return hit_phi;
  }

  inline void setGeomParams(int chamberID, int layerID, double layerphi, double R, double eps, double offset) {

    layer_params.insert(std::pair<vID,LAYER>(vID(chamberID,layerID),LAYER(layerphi,R,eps,offset)));

   }

  inline void setWiresInLayer(int chamber, int layer, int numWires)
  {
    updateParams(chamber,layer);
    for (int i = 0; i<numWires; ++i) {

      auto phi_start = _currentLayerphi * (i+0.5) + m_offset;
      auto phi_end = phi_start + _currentLayerphi;

      TVector3 Wstart = returnWirePosition(phi_start, 1);
      TVector3 Wend = returnWirePosition(phi_end, -1);

      TVector3 Wmid = (Wstart+Wend)*(1/2.0);
      TVector3 Wdirection = (Wend - Wstart);

      m_wiresPositions[layer].push_back(std::make_pair(Wmid, Wdirection));
      }
  }

  inline auto returnAllWires() const { return m_wiresPositions; }

//  inline TVector3 returnWirePosition(double angle, int sign) const {
//    TVector3 w(0, 0, 0);
//    w.SetX(_currentRadius * std::cos(angle));
//    w.SetY(_currentRadius * std::sin(angle));
//    w.SetZ(sign * m_detectorLength / 2.0);
//    return w;
//  }

  void updateParams(int chamber, int layer)  const{
    auto it_end = layer_params.cend();
    --it_end;
    double LayerPhi = it_end->second.layerphi;
    double Radius = it_end->second.R;
    double Eps = it_end->second.eps;
    double Offset = it_end->second.offset;

    CID v1(chamber,layer);
    auto map_it = layer_params.find(v1);
    if (map_it != layer_params.cend()) {
     LayerPhi = map_it->second.layerphi;
     Radius = map_it->second.R;
     Eps = map_it->second.eps;
     Offset = map_it->second.offset;
    } else { std::cout << " Sorry, pair with  key " << layer << " not in map" << std::endl; }

    _currentLayerphi = LayerPhi;
    _currentRadius = Radius;
    m_epsilon = Eps;
    m_offset = Offset;
 }

protected:

  double phi(const CellID& cID) const;
  std::map<vID,LAYER> layer_params; // <{chamberID,layerID}, {layerphi, R, eps, offset}>
  std::map<int, std::vector<std::pair<TVector3, TVector3> >> m_wiresPositions; // < layer, vec<WireMidpoint, WireDirection> >

  inline TVector3 returnWirePosition(double angle, int sign) const {
    TVector3 w(0, 0, 0);
    w.SetX(_currentRadius * std::cos(angle));
    w.SetY(_currentRadius * std::sin(angle));
    w.SetZ(sign * m_detectorLength / 2.0);
    return w;
  }

  inline double returnAlpha() const {
    double alpha = 2 * std::asin(m_detectorLength * std::tan(m_epsilon0)/(2 * _currentRadius));
    return alpha;
 }

  double m_cellSize;
  double m_epsilon0;
  double m_detectorLength;
  double m_DC_inner_rbegin;
  double m_DC_inner_rend;
  double m_DC_outer_rbegin;
  double m_DC_outer_rend;
  double m_DC_inner_rmin;
  double m_DC_inner_rmax;
  double m_DC_outer_rmin;
  double m_DC_outer_rmax;
  double m_layer_width;
  double m_safe_distance;
  int m_DC_inner_layer_number;
  int m_DC_outer_layer_number;

  std::string m_phiID;
  std::string layer_id;

  // Current parameters of the layer: sizePhi
  mutable double _currentLayerphi;
  mutable double _currentRadius;
  mutable double m_epsilon;
  mutable double m_offset;

};
}
}
#endif /* DETSEGMENTATION_GRIDDRIFTCHAMBER_H */
