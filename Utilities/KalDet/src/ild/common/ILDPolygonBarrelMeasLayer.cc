
#include "ILDPolygonBarrelMeasLayer.h"

#include <UTIL/BitField64.h>
#include <UTIL/ILDConf.h>

#include "kaltest/TVTrack.h"
#include "TVector3.h"
#include "TMath.h"
#include "TRotMatrix.h"
#include "TBRIK.h"
#include "TNode.h"
#include "TString.h"

//#include <EVENT/TrackerHitPlane.h>
#include <edm4hep/TrackerHit.h>

#include <math.h>
#include <assert.h>
#include <algorithm>

// #include "streamlog/streamlog.h"


ILDPolygonBarrelMeasLayer::ILDPolygonBarrelMeasLayer(TMaterial &min,
                                                     TMaterial &mout,
                                                     double   Bz,
                                                     double   sortingPolicy,
                                                     double   r0,       // min distance to the z-axis
                                                     double   lhalf,    // half length
                                                     int      nsides,   
                                                     double   zpos,     // z of the centre 
                                                     double   phi0,     // phi of the first normal following the xaxis positive rotation
                                                     std::vector<int>      CellIDs,                                                     
                                                     bool     is_active,
                                                     const Char_t    *name) : 
ILDVMeasLayer(min, mout, Bz, CellIDs, is_active, name), _sortingPolicy(0.),
_r0(r0),_lhalf(lhalf),_nsides(nsides),_zpos(zpos),_phi0(phi0)

{
    
  _segment_dphi = 2.0*M_PI / _nsides; 
  
  phi0 = angular_range_2PI(phi0);
  
  _start_phi = phi0 - 0.5*_segment_dphi;
 
  _start_phi = angular_range_2PI(_start_phi);

  for (int i = 0; i < _nsides; ++i) {
    
    const TVector3  center(r0*cos(phi0),r0*sin(phi0),zpos);
    const TVector3  normal(center);
    double width = 2*(_r0*tan(_segment_dphi*0.5));
    
    double phi = i*_segment_dphi+_start_phi;
    
    ILDParallelPlanarMeasLayer p(min,mout,r0,phi,Bz,sortingPolicy,width,2*lhalf,0.0,0.0,0.0,false); // NOTE: No Physical Offset or Coordinate Offset allowed
    _planes.push_back(p);
    
  }

  
  _rmax = r0*sin((_r0*tan(_segment_dphi*0.5)));
  
  _enclosing_cylinder = new  TCylinder(_rmax, lhalf, 0., 0., 0.);
  
  
}


TKalMatrix ILDPolygonBarrelMeasLayer::XvToMv(const TVector3 &xv) const
{
  
  // Calculate measurement vector (hit coordinates) from global coordinates:
  int index = this->get_plane_index(xv.Phi());
  
  return _planes[index].XvToMv(xv);
  
}


TVector3 ILDPolygonBarrelMeasLayer::HitToXv(const TVTrackHit &vht) const
{
  
  //SJA:FIXME: in order to use real local coordinates we would have to get the CELLID from the ILDPlanarHit, this would tell us in which segment the hit was in 
  
  // streamlog_out(ERROR) << "ILDPolygonBarrelMeasLayer::HitToXv Not implemented: exit(1) called from " << __FILE__ << "   line " << __LINE__ << std::endl; 
  exit(1);

  return TVector3(0.,0.,0.);
  
}

void ILDPolygonBarrelMeasLayer::CalcDhDa(const TVTrackHit &vht,
                                         const TVector3   &xxv,
                                         const TKalMatrix &dxphiada,
                                         TKalMatrix &H)  const
{

    //SJA:FIXME: in order to use real local coordinates we would have to get the CELLID from the ILDPlanarHit, this would tell us in which segment the hit was in 
  
  // streamlog_out(ERROR) << "ILDPolygonBarrelMeasLayer::CalcDhDa Not implemented: exit(1) called from " << __FILE__ << "   line " << __LINE__ << std::endl; 
  exit(1);
  
}


Int_t ILDPolygonBarrelMeasLayer::CalcXingPointWith(const TVTrack  &hel,
                                                    TVector3 &xx,
                                                    Double_t &phi,
                                                    Int_t     mode,
                                                    Double_t  eps) const{

  int crosses = _enclosing_cylinder->CalcXingPointWith(hel, xx, phi, mode);
  
  if (crosses == 0) {
    return 0;
  }
  
  int index = this->get_plane_index(xx.Phi());
  
  return _planes[index].CalcXingPointWith(hel,xx,phi,mode,eps);

    
}

Bool_t ILDPolygonBarrelMeasLayer::IsOnSurface(const TVector3 &xx) const
{
  
  //  streamlog_out(DEBUG0) << "IsOnSurface " << std::endl;  
  int index = this->get_plane_index(xx.Phi());
  
  return _planes[index].IsOnSurface(xx); 
    
}


ILDVTrackHit* ILDPolygonBarrelMeasLayer::ConvertLCIOTrkHit(edm4hep::ConstTrackerHit trkhit) const {
  
  // streamlog_out(ERROR) << "ILDPolygonBarrelMeasLayer::ConvertLCIOTrkHit Not implemented: exit(1) called from " << __FILE__ << "   line " << __LINE__ << std::endl; 
  exit(1);
    
}

bool ILDPolygonBarrelMeasLayer::IsOutside(const TVector3 &xx) const
{
  int index = this->get_plane_index(xx.Phi());
  return _planes[index].IsOutside(xx); 
}

double ILDPolygonBarrelMeasLayer::CalcS(const TVector3 &xx) const {

  int index = this->get_plane_index(xx.Phi());
  return _planes[index].CalcS(xx); 
  
}

TMatrixD ILDPolygonBarrelMeasLayer::CalcDSDx(const TVector3 &xx) const {
  
  int index = this->get_plane_index(xx.Phi());
  return _planes[index].CalcDSDx(xx); 
  
}

/** Get the intersection and the CellID, needed for multilayers */
int ILDPolygonBarrelMeasLayer::getIntersectionAndCellID(const TVTrack  &hel,
                                                        TVector3 &xx,
                                                        Double_t &phi,
                                                        Int_t    &CellID,
                                                        Int_t     mode,
                                                        Double_t  eps) const {
  
  
  int crosses = this->CalcXingPointWith(hel, xx, phi, mode, eps);
  
  if ( crosses != 0 ) {
    
    unsigned int plane = this->get_plane_index( xx.Phi() );
    
    lcio::BitField64 bf(  UTIL::ILDCellID0::encoder_string ) ;
    bf.setValue( this->getCellIDs()[0] ) ; // get the first cell_ID, this will have the correct sensor value
    
    bf[lcio::ILDCellID0::module] = plane;
    CellID = bf.lowWord();
    
  }
  
  return crosses;
  
}


unsigned int ILDPolygonBarrelMeasLayer::get_plane_index(double phi) const {

  phi = angular_range_2PI(phi-_start_phi);
  return unsigned(floor(phi/_segment_dphi));
  
}


double ILDPolygonBarrelMeasLayer::angular_range_2PI( double phi ) const {
  
  //bring phi_point into range 0 < phi < +2PI
  while (phi < 0) {
    phi += 2.0 * M_PI;
  }
  while (phi >= 2.0*M_PI) {
    phi -= 2.0 * M_PI;
  }
  
  return phi;
  
}

