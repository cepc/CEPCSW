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
  registerIdentifier("layerID", "layer id", layer_id, "layer");
}

GridDriftChamber::GridDriftChamber(const BitFieldCoder* decoder) : Segmentation(decoder) {

  _type = "GridDriftChamber";
  _description = "Drift chamber segmentation in the global coordinates";

  registerParameter("cell_size", "cell size", m_cellSize, 1., SegmentationParameter::LengthUnit);
  registerParameter("detector_length", "Length of the wire", m_detectorLength, 1., SegmentationParameter::LengthUnit);
  registerIdentifier("identifier_phi", "Cell ID identifier for phi", m_phiID, "cellID");
  registerIdentifier("layerID", "layer id", layer_id, "layer");
  registerParameter("layer_width", "layer_width", m_layer_width, 0., SegmentationParameter::LengthUnit);
  registerParameter("DC_rbegin", "DC_rbegin", m_DC_rbegin, 0., SegmentationParameter::LengthUnit);
  registerParameter("DC_rend", "DC_rend", m_DC_rend, 0., SegmentationParameter::LengthUnit);
}

Vector3D GridDriftChamber::position(const CellID& /*cID*/) const {
  Vector3D cellPosition = {0, 0, 0};
  return cellPosition;
}

CellID GridDriftChamber::cellID(const Vector3D& /*localPosition*/, const Vector3D& globalPosition,
                                const VolumeID& vID) const {

  CellID cID = vID;

  int chamberID = _decoder->get(cID, "chamber");
  int layerid = _decoder->get(cID, "layer");

  double posx = globalPosition.X;
  double posy = globalPosition.Y;
  double posz = globalPosition.Z;

  updateParams(chamberID,layerid);

  TVector3 Phi0 = returnPhi0(chamberID,layerid,posz);
  double phi0 = phiFromXY2(Phi0);

  double phi_hit = phiFromXY(globalPosition);
  double offsetphi= m_offset;
  double _lphi;

  _lphi = phi_hit - phi0 + _currentLayerphi/2.;
  if(_lphi<0.){
      _lphi+=2*M_PI;
  }else if(_lphi>2*M_PI){
      _lphi=fmod(_lphi,2*M_PI);
  }
  int cellID=floor(_lphi/_currentLayerphi);

  _decoder->set(cID, m_phiID, cellID);

  return cID;
}

double GridDriftChamber::phi(const CellID& cID) const {
    CellID phiValue = _decoder->get(cID, m_phiID);
    return binToPosition(phiValue, _currentLayerphi, m_offset);
}

void GridDriftChamber::cellposition(const CellID& cID, TVector3& Wstart,
        TVector3& Wend) const {

    auto chamberIndex = _decoder->get(cID, "chamber");
    auto layerIndex = _decoder->get(cID, "layer");
    updateParams(chamberIndex,layerIndex);

    double phi_start = phi(cID);
    double phi_end = phi_start + returnAlpha();

    Wstart = returnWirePosition(phi_start, -1);
    Wend = returnWirePosition(phi_end, 1);
}

TVector3 GridDriftChamber::returnPhi0(int chamber,int layer, double z) const
{
    updateParams(chamber,layer);

    double phi_wire_start = binToPosition(0 , _currentLayerphi, m_offset);
    double phi_wire_end = phi_wire_start + returnAlpha();

    TVector3 wire_start = returnWirePosition(phi_wire_start, -1);
    TVector3 wire_end = returnWirePosition(phi_wire_end, 1);

    double ratio = fabs(z - wire_start.Z())/fabs(wire_end.Z() - wire_start.Z());
    double x_pos = ratio * (wire_end.X() - wire_start.X()) + wire_start.X();
    double y_pos = ratio * (wire_end.Y() - wire_start.Y()) + wire_start.Y();

    return TVector3(x_pos,y_pos,z);
}


void GridDriftChamber::cellposition2(int chamber,int layer, int cell,
        TVector3& Wstart, TVector3& Wend) const {
    updateParams(chamber,layer);
    double phi_start = binToPosition(cell, _currentLayerphi, m_offset);
    double phi_end = phi_start + returnAlpha();

    Wstart = returnWirePosition(phi_start, -1);
    Wend = returnWirePosition(phi_end, 1);
}

TVector3 GridDriftChamber::LineLineIntersect(TVector3& p1, TVector3& p2, TVector3& p3, TVector3& p4) const {

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

double GridDriftChamber::distanceTrackWire(const CellID& cID, const TVector3& hit_start,
        const TVector3& hit_end) const {

    TVector3 Wstart = {0,0,0};
    TVector3 Wend = {0,0,0};
    cellposition(cID,Wstart,Wend);

    TVector3 a = (hit_end - hit_start).Unit();
    TVector3 b = (Wend - Wstart).Unit();
    TVector3 c = Wstart - hit_start;

    double num = std::abs(c.Dot(a.Cross(b)));
    double denum = (a.Cross(b)).Mag();

    double DCA = 0;

    if (denum) {
        DCA = num / denum;
    }

    return DCA;
}

double GridDriftChamber::distanceTrackWire2(const CellID& cID, const TVector3& hit_pos) const {

    TVector3 Wstart = {0,0,0};
    TVector3 Wend = {0,0,0};
    cellposition(cID,Wstart,Wend);

    TVector3 denominator = Wend - Wstart;
    TVector3  numerator = denominator.Cross(Wstart-hit_pos);

    double DCA = numerator.Mag()/denominator.Mag() ;

    return DCA;
}

TVector3 GridDriftChamber::distanceClosestApproach(const CellID& cID, const TVector3& hitPos, TVector3& PCA) const {
    // Distance of the closest approach between a single hit (point) and the closest wire

    TVector3 Wstart = {0,0,0};
    TVector3 Wend = {0,0,0};
    cellposition(cID,Wstart,Wend);

    TVector3 temp = (Wend + Wstart);
    TVector3 Wmid(temp.X() / 2.0, temp.Y() / 2.0, temp.Z() / 2.0);

    double hitPhi = hitPos.Phi();
    if (hitPhi < 0) {
        hitPhi = hitPhi + 2 * M_PI;
    }

    PCA = Wstart + ((Wend - Wstart).Unit()).Dot((hitPos - Wstart)) * ((Wend - Wstart).Unit());
    TVector3 dca = hitPos - PCA;

    return dca;
}

TVector3 GridDriftChamber::Line_TrackWire(const CellID& cID, const TVector3& hit_start, const TVector3& hit_end) const {
    // The line connecting a particle track to the closest wire
    // Returns the vector connecting the both
    TVector3 Wstart = {0,0,0};
    TVector3 Wend = {0,0,0};
    cellposition(cID,Wstart,Wend);

    TVector3 P1 = hit_start;
    TVector3 P2 = hit_end;
    TVector3 P3 = Wstart;
    TVector3 P4 = Wend;

    TVector3 intersect = LineLineIntersect(P1, P2, P3, P4);
    return intersect;
}

// Get the wire position for a z
TVector3 GridDriftChamber::wirePos_vs_z(const CellID& cID, const double& zpos) const {

    TVector3 Wstart = {0,0,0};
    TVector3 Wend = {0,0,0};
    cellposition(cID,Wstart,Wend);

    double t = (zpos - Wstart.Z())/(Wend.Z()-Wstart.Z());
    double x = Wstart.X()+t*(Wend.X()-Wstart.X());
    double y = Wstart.Y()+t*(Wend.Y()-Wstart.Y());

    TVector3 wireCoord(x, y, zpos);
    return wireCoord;
}

TVector3 GridDriftChamber::IntersectionTrackWire(const CellID& cID, const TVector3& hit_start, const TVector3& hit_end) const {
    // Intersection between the particle track and the wire assuming that the track between hit_start and hit_end is linear

    TVector3 Wstart = {0,0,0};
    TVector3 Wend = {0,0,0};
    cellposition(cID,Wstart,Wend);

    TVector3 P1 = hit_start;
    TVector3 V1 = hit_end-hit_start;

    TVector3 P2 = Wstart;
    TVector3 V2 = Wend - Wstart;

    TVector3 denom = V1.Cross(V2);
    double mag_denom = denom.Mag();

    TVector3 intersect(0, 0, 0);

    if (mag_denom !=0)
    {
        TVector3 num = ((P2-P1)).Cross(V2);
        double mag_num = num.Mag();
        double a = mag_num / mag_denom;

        intersect = P1 + a * V1;

    }
    return intersect;
}

double GridDriftChamber::Distance(const CellID& cID, const TVector3& pointIn, const TVector3& pointOut, TVector3& hitPosition, TVector3& PCA) const {

    //For two lines r=r1+t1.v1 & r=r2+t2.v2
    //the closest approach is d=|(r2-r1).(v1 x v2)|/|v1 x v2|
    //the point where closest approach are
    //t1=(v1 x v2).[(r2-r1) x v2]/[(v1 x v2).(v1 x v2)]
    //t2=(v1 x v2).[(r2-r1) x v1]/[(v1 x v2).(v1 x v2)]
    //if v1 x v2=0 means two lines are parallel
    //d=|(r2-r1) x v1|/|v1|

    double t1,distance,dInOut,dHitIn,dHitOut;
    //Get wirepoint @ endplate
    TVector3 west = {0,0,0};
    TVector3 east = {0,0,0};
    cellposition(cID,west,east);
    TVector3 wireLine=east - west;
    TVector3 hitLine=pointOut - pointIn;

    TVector3 hitXwire=hitLine.Cross(wireLine);
    TVector3 wire2hit=east-pointOut;
    //Hitposition is the position on hit line where closest approach
    //of two lines, but it may out the area from pointIn to pointOut
    if(hitXwire.Mag()==0){
        distance=wireLine.Cross(wire2hit).Mag()/wireLine.Mag();
        hitPosition=pointIn;
    }else{
        t1=hitXwire.Dot(wire2hit.Cross(wireLine))/hitXwire.Mag2();
        hitPosition=pointOut+t1*hitLine;

        dInOut=(pointOut-pointIn).Mag();
        dHitIn=(hitPosition-pointIn).Mag();
        dHitOut=(hitPosition-pointOut).Mag();
        if(dHitIn<=dInOut && dHitOut<=dInOut){ //Between point in & out
            distance=fabs(wire2hit.Dot(hitXwire)/hitXwire.Mag());
        }else if(dHitOut>dHitIn){ // out pointIn
            distance=wireLine.Cross(pointIn-east).Mag()/wireLine.Mag();
            hitPosition=pointIn;
        }else{ // out pointOut
            distance=wireLine.Cross(pointOut-east).Mag()/wireLine.Mag();
            hitPosition=pointOut;
        }
    }

    PCA = west + ((east - west).Unit()).Dot((hitPosition - west)) * ((east - west).Unit());

    return distance;
}

}
}
