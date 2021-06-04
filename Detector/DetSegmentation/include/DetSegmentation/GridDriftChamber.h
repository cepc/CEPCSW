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
  virtual TVector3 distanceClosestApproach(const CellID& cID, const TVector3& hitPos) const;
  virtual TVector3 Line_TrackWire(const CellID& cID, const TVector3& hit_start, const TVector3& hit_end) const;
  virtual TVector3 IntersectionTrackWire(const CellID& cID, const TVector3& hit_start, const TVector3& hit_end) const;
  virtual TVector3 wirePos_vs_z(const CellID& cID, const double& zpos) const;

//  double phi(const CellID& cID) const;
  inline double cell_Size() const { return m_cellSize; }
  inline double epsilon0() const { return m_epsilon0; }
  inline double detectorLength() const { return m_detectorLength; }
  inline double safe_distance() const { return m_safe_distance; }
  inline double layer_width() const { return m_layer_width; }
  inline double DC_rbegin() const { return m_DC_rbegin; }
  inline double DC_rend() const { return m_DC_rend; }
  inline double DC_rmax() const { return m_DC_rmax; }
  inline double DC_rmin() const { return m_DC_rmin; }
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

  TVector3 LineLineIntersect(TVector3 p1, TVector3 p2, TVector3 p3, TVector3 p4) const {
    TVector3 p13, p43, p21;
    double d1343, d4321, d1321, d4343, d2121;
    double numer, denom;
    double mua, mub;
    TVector3 pa, pb;

    p13.SetX(p1.X() - p3.X());
    p13.SetY(p1.Y() - p3.Y());
    p13.SetZ(p1.Z() - p3.Z());
    p43.SetX(p4.X() - p3.X());
    p43.SetY(p4.Y() - p3.Y());
    p43.SetZ(p4.Z() - p3.Z());
    /* if (ABS(p43.X())  < EPS && ABS(p43.Y())  < EPS && ABS(p43.Z())  < EPS) */
    /*    return(FALSE); */
    p21.SetX(p2.X() - p1.X());
    p21.SetY(p2.Y() - p1.Y());
    p21.SetZ(p2.Z() - p1.Z());
    /* if (ABS(p21.X())  < EPS && ABS(p21.Y())  < EPS && ABS(p21.Z())  < EPS) */
    /*    return(FALSE); */

    d1343 = p13.X() * p43.X() + p13.Y() * p43.Y() + p13.Z() * p43.Z();
    d4321 = p43.X() * p21.X() + p43.Y() * p21.Y() + p43.Z() * p21.Z();
    d1321 = p13.X() * p21.X() + p13.Y() * p21.Y() + p13.Z() * p21.Z();
    d4343 = p43.X() * p43.X() + p43.Y() * p43.Y() + p43.Z() * p43.Z();
    d2121 = p21.X() * p21.X() + p21.Y() * p21.Y() + p21.Z() * p21.Z();

    denom = d2121 * d4343 - d4321 * d4321;
    /* if (ABS(denom) < EPS) */
    /*    return(FALSE); */
    numer = d1343 * d4321 - d1321 * d4343;

    mua = numer / denom;
    mub = (d1343 + d4321 * (mua)) / d4343;

    pa.SetX(p1.X() + mua * p21.X());
    pa.SetY(p1.Y() + mua * p21.Y());
    pa.SetZ(p1.Z() + mua * p21.Z());
    pb.SetX(p3.X() + mub * p43.X());
    pb.SetY(p3.Y() + mub * p43.Y());
    pb.SetZ(p3.Z() + mub * p43.Z());

    return pb - pa;
  }


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
  double m_layer_width;
  double m_safe_distance;
  double m_DC_rbegin;
  double m_DC_rend;
  double m_DC_rmax;
  double m_DC_rmin;

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
