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
  virtual double distanceTrackWire2(const CellID& cID, const TVector3& hit_pos) const;
  virtual void cellposition(const CellID& cID, TVector3& Wstart, TVector3& Wend) const;
  virtual void cellposition2(int chamber, int layer, int cell, TVector3& Wstart, TVector3& Wend) const;
  TVector3 LineLineIntersect(TVector3& p1, TVector3& p2, TVector3& p3, TVector3& p4) const;
  virtual TVector3 distanceClosestApproach(const CellID& cID, const TVector3& hitPos, TVector3& PCA) const;
  virtual TVector3 Line_TrackWire(const CellID& cID, const TVector3& hit_start, const TVector3& hit_end) const;
  virtual TVector3 IntersectionTrackWire(const CellID& cID, const TVector3& hit_start, const TVector3& hit_end) const;
  virtual TVector3 wirePos_vs_z(const CellID& cID, const double& zpos) const;
  virtual double Distance(const CellID& cID, const TVector3& pointIn, const TVector3& pointOut, TVector3& hitPosition, TVector3& PCA) const;
  virtual TVector3 returnPhi0(int chamber,int layer, double z) const;


//  double phi(const CellID& cID) const;
  inline double cell_Size() const { return m_cellSize; }
  inline double epsilon() const { return m_epsilon; }
  inline double detectorLength() const { return m_detectorLength; }
  inline double layer_width() const { return m_layer_width; }
  inline double DC_rbegin() const { return m_DC_rbegin; }
  inline double DC_rend() const { return m_DC_rend; }
  inline const std::string& fieldNamePhi() const { return m_phiID; }
  inline const std::string& Layerid() const { return layer_id; }

  // Setters

  inline double phiFromXY(const Vector3D& aposition) const {
    double hit_phi =  std::atan2(aposition.Y, aposition.X) ;
    if( hit_phi < 0 ) { hit_phi += 2 * M_PI; }
    return hit_phi;
  }

 inline double phiFromXY2(const TVector3& aposition) const {
    double hit_phi =  std::atan2(aposition.Y(), aposition.X()) ;
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

  inline Vector3D returnPosWire0(double z) const {
    double alpha = returnAlpha();
    double t = 0.5 * (1 - 2.0 * z / m_detectorLength);
    double x = _currentRadius * (1 + t * (std::cos(alpha) - 1));
    double y = _currentRadius * t * std::sin(alpha);

    Vector3D vec(x, y, z);
    return vec;
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
    double alpha = 2 * std::asin(m_detectorLength * std::tan(m_epsilon)/(2 * _currentRadius));
    return alpha;
 }

  double m_cellSize;
  double m_detectorLength;
  double m_layer_width;
  double m_DC_rbegin;
  double m_DC_rend;

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
