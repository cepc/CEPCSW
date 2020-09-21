//*************************************************************************
//* ===================
//*  ILDPlanarMeasLayer Class
//* ===================
//*
//* (Description)
//*   Sample measurement layer class used by ILDPlanarHit.
//* (Requires)
//*     ILDVMeasLayer
//* (Provides)
//*     class ILDPlanarMeasLayer
//*
//*************************************************************************
//
#include <iostream>
#include <cmath>

#include "ILDPlanarMeasLayer.h"
#include "ILDPlanarHit.h"

#include "kaltest/TVTrack.h"
#include "TVector3.h"
#include "TMath.h"
#include "TRotMatrix.h"
#include "TBRIK.h"
#include "TNode.h"
#include "TString.h"

#include "gearimpl/Vector3D.h"

//#include <EVENT/TrackerHitPlane.h>

// #include "streamlog/streamlog.h"

ILDPlanarMeasLayer::ILDPlanarMeasLayer(TMaterial &min,
                                       TMaterial &mout,
                                       const TVector3  &center,
                                       const TVector3  &normal,
                                       Double_t   Bz,
                                       Double_t   SortingPolicy,
                                       Double_t   xiwidth,
                                       Double_t   zetawidth,
                                       Double_t   xioffset,
                                       Double_t   UOrigin,
                                       Bool_t     is_active,
                                       Int_t      CellID,
                                       const Char_t    *name)
: ILDVMeasLayer(min, mout, Bz, is_active, CellID, name),
TPlane(center, normal),
fSortingPolicy(SortingPolicy),
fXiwidth(xiwidth),
fZetawidth(zetawidth),
fXioffset(xioffset),
fUOrigin(UOrigin)

{
  
  // streamlog_out(DEBUG0) << "ILDPlanarMeasLayer created" 
  // << " Layer x0 = " << this->GetXc().X() 
  // << " y0 = " << this->GetXc().Y() 
  // << " z0 = " << this->GetXc().Z() 
  // << " R = " << this->GetXc().Perp() 
  // << " phi = " << this->GetXc().Phi() 
  // << " xiwidth = " << fXiwidth 
  // << " zetawidth = " << fZetawidth 
  // << " xioffset = " << fXioffset 
  // << " UOrigin = " << UOrigin
  // << " is_active = " << is_active 
  // << " CellID = " << CellID 
  // << " name = " << this->ILDVMeasLayer::GetName()  
  // << std::endl ;
  
  
}

ILDPlanarMeasLayer::~ILDPlanarMeasLayer()
{
}

TKalMatrix ILDPlanarMeasLayer::XvToMv(const TVector3 &xv) const
{
  // Calculate hit coordinate information:
  //    mv(0,0) = xi 
  //     (1,0) = zeta
  
//  streamlog_out(DEBUG0) << "\t ILDPlanarMeasLayer::XvToMv: "
//  << " x = " << xv.X() 
//  << " y = " << xv.Y() 
//  << " z = " << xv.Z() 
//  << std::endl;
  
  TKalMatrix mv(kMdim,1);
  
  double cos_phi = GetNormal().X()/GetNormal().Perp();
  double sin_phi = GetNormal().Y()/GetNormal().Perp();
  
  double delta_x = xv.X() - GetXc().X();
  double delta_y = xv.Y() - GetXc().Y();
  double delta_z = xv.Z() - GetXc().Z();

  double delta_t = (delta_x * sin_phi - delta_y * cos_phi) ;
  
  mv(0,0) = delta_t + fUOrigin;

  mv(1,0) = delta_z ;
  
//  streamlog_out(DEBUG0) << "\t ILDPlanarMeasLayer::XvToMv: "
//  << " mv(0,0) = " << mv(0,0)
//  << " mv(1,0) = " << mv(1,0)
//  << std::endl;

  
  return mv;
}

TKalMatrix ILDPlanarMeasLayer::XvToMv(const TVTrackHit &,
                                      const TVector3   &xv) const
{
  return XvToMv(xv);
}

TVector3 ILDPlanarMeasLayer::HitToXv(const TVTrackHit &vht) const
{
 
  
//  streamlog_out(DEBUG0) << "\t ILDPlanarMeasLayer::HitToXv: "
//  << " vht(0,0) = " << vht(0,0)
//  << " vht(1,0) = " << vht(1,0)
//  << std::endl;
  

  //const ILDPlanarHit &ht = dynamic_cast<const ILDPlanarHit &>(vht);
  
  double cos_phi = GetNormal().X()/GetNormal().Perp();
  double sin_phi = GetNormal().Y()/GetNormal().Perp();

  double delta_t = vht(0,0) - fUOrigin;
  
  double x =  delta_t * sin_phi + this->GetXc().X();
  double y = -delta_t * cos_phi + this->GetXc().Y();
  
  double z =  vht(1,0) + this->GetXc().Z();
  
//  streamlog_out(DEBUG0) << "\t ILDPlanarMeasLayer::HitToXv: "
//  << " x = " << x 
//  << " y = " << y 
//  << " z = " << z 
//  << std::endl;

  
  return TVector3(x,y,z);
}

void ILDPlanarMeasLayer::CalcDhDa(const TVTrackHit &vht,
                                  const TVector3   &xxv,
                                  const TKalMatrix &dxphiada,
                                  TKalMatrix &H)  const
{
  // Calculate
  //    H = (@h/@a) = (@phi/@a, @z/@a)^t
  // where
  //        h(a) = (phi, z)^t: expected meas vector
  //        a = (drho, phi0, kappa, dz, tanl, t0)
  //
  
  double cos_phi = GetNormal().X()/GetNormal().Perp();
  double sin_phi = GetNormal().Y()/GetNormal().Perp();
  
  
  Int_t sdim = H.GetNcols();
  Int_t hdim = TMath::Max(5,sdim-1);
  
  // Set H = (@h/@a) = (@d/@a, @z/@a)^t
  
  for (Int_t i=0; i<hdim; i++) {

    H(0,i) =  sin_phi * dxphiada(0,i) - cos_phi * dxphiada(1,i) ;   
    H(1,i) =  dxphiada(2,i);
    
  }
  if (sdim == 6) {
    H(0,sdim-1) = 0.;
    H(1,sdim-1) = 0.;
  }
  
}


//#define DEBUG_ISONSURFACE 1 
Bool_t ILDPlanarMeasLayer::IsOnSurface(const TVector3 &xx) const
{
  Double_t xi   = (xx.Y()-GetXc().Y())*GetNormal().X()/GetNormal().Perp() - (xx.X()-GetXc().X())*GetNormal().Y()/GetNormal().Perp() ;
  Double_t zeta = xx.Z() - GetXc().Z();
    
  bool onSurface = false ;
  
  if( (xx.X()-GetXc().X())*GetNormal().X() + (xx.Y()-GetXc().Y())*GetNormal().Y() < 1e-4){
    if( xi <= GetXioffset() + GetXiwidth()/2  && xi >= GetXioffset() - GetXiwidth()/2  && TMath::Abs(zeta) <= GetZetawidth()/2){
      onSurface = true;
    }
    
#ifdef DEBUG_ISONSURFACE
    else{
      // streamlog_out(DEBUG4) << "ILDPlanarMeasLayer::IsOnSurface: Point not within boundary x = " << xx.x() << " y = " << xx.y() << " z = " << xx.z() << " r = " << xx.Perp() << " phi = " << xx.Phi() << std::endl;   
      // streamlog_out(DEBUG4) << "xi = " << xi << " xi_max = " << GetXioffset() + GetXiwidth()/2 << " xi_min = " << GetXioffset() - GetXiwidth()/2 << " zeta = " << zeta << " zeta_min = " << -GetZetawidth()/2 << " zeta_max " << GetZetawidth()/2 << " Xioffset = " << GetXioffset() << std::endl;     
      onSurface = false;
      // streamlog_out(DEBUG4) << " xi <= GetXioffset() + GetXiwidth()/2 = " << (xi <= GetXioffset() + GetXiwidth()/2) << std::endl ;
      // streamlog_out(DEBUG4) << " xi >= GetXioffset() - GetXiwidth()/2 = " << (xi >= GetXioffset() - GetXiwidth()/2) << std::endl;
      // streamlog_out(DEBUG4) << " TMath::Abs(zeta) <= GetZetawidth()/2 = " << (TMath::Abs(zeta) <= GetZetawidth()/2) << std::endl;
    }
#endif

  }

#ifdef DEBUG_ISONSURFACE  
  else{
    // streamlog_out(DEBUG4) << "ILDPlanarMeasLayer::IsOnSurface: Point not on surface x = " << xx.x() << " y = " << xx.y() << " z = " << xx.z() << " r = " << xx.Perp() << " phi = " << xx.Phi() << std::endl;   
    // streamlog_out(DEBUG4) << "Distance from plane " << (xx.X()-GetXc().X())*GetNormal().X() + (xx.Y()-GetXc().Y())*GetNormal().Y() << std::endl;
  }
  if( onSurface == false ) {
    // streamlog_out(DEBUG) << "x0 " <<  GetXc().X() << std::endl;  
    // streamlog_out(DEBUG) << "y0 " <<  GetXc().Y() << std::endl;  
    // streamlog_out(DEBUG) << "z0 " <<  GetXc().Z() << std::endl;  
    // streamlog_out(DEBUG) << "GetNormal().X() " <<  GetNormal().X() << std::endl;  
    // streamlog_out(DEBUG4) << "GetNormal().Y() " <<  GetNormal().Y() << std::endl;  
    // streamlog_out(DEBUG4) << "GetNormal().Perp() " << GetNormal().Perp() << std::endl;  
    // streamlog_out(DEBUG4) << "GetNormal().X()/GetNormal().Perp() " << GetNormal().X()/GetNormal().Perp() << std::endl;  
    // streamlog_out(DEBUG4) << "GetNormal().Y()/GetNormal().Perp() " << GetNormal().Y()/GetNormal().Perp() << std::endl;  
    // streamlog_out(DEBUG4) << "xx.X()-GetXc().X() " <<  xx.X()-GetXc().X() << std::endl;  
    // streamlog_out(DEBUG4) << "xx.Y()-GetXc().Y() " <<  xx.Y()-GetXc().Y() << std::endl;  
    
    // streamlog_out(DEBUG4) << "zeta " << zeta << std::endl;  
    // streamlog_out(DEBUG4) << "xi "   << xi   << std::endl;  
    // streamlog_out(DEBUG4) << "zeta half width " << GetZetawidth()/2 << std::endl;  
    // streamlog_out(DEBUG4) << "xi half width " << GetXiwidth()/2 << std::endl;  
    // streamlog_out(DEBUG4) << "offset  " << GetXioffset() << std::endl;  
    
    // streamlog_out(DEBUG4) << "distance from plane " << (xx.X()-GetXc().X())*GetNormal().X() + (xx.Y()-GetXc().Y())*GetNormal().Y() << std::endl; 
  }
#endif
  
  return onSurface;
  
}


ILDVTrackHit* ILDPlanarMeasLayer::ConvertLCIOTrkHit(edm4hep::ConstTrackerHit trkhit) const {
  //std::cout << "ILDPlanarMeasLayer::ConvertLCIOTrkHit " << trkhit << " type=" << trkhit.getType() << std::endl;
  //EVENT::TrackerHitPlane* plane_hit = dynamic_cast<EVENT::TrackerHitPlane*>( trkhit ) ;
  if((trkhit.getType()&8)!=8){
    std::cout << "ILDPlanarMeasLayer::ConvertLCIOTrkHit Warning: type is not 8, but " << (trkhit.getType()&8) << std::endl;
    return NULL;
  }
  //if( plane_hit == NULL )  return NULL; // SJA:FIXME: should be replaced with an exception  

  //gear::Vector3D U(1.0,plane_hit.getU()[1],plane_hit.getU()[0],gear::Vector3D::spherical);
  //gear::Vector3D V(1.0,plane_hit.getV()[1],plane_hit.getV()[0],gear::Vector3D::spherical);
  gear::Vector3D U(1.0,trkhit.getCovMatrix(1),trkhit.getCovMatrix(0),gear::Vector3D::spherical);
  gear::Vector3D V(1.0,trkhit.getCovMatrix(4),trkhit.getCovMatrix(3),gear::Vector3D::spherical);
  gear::Vector3D Z(0.0,0.0,1.0);
  
  const float eps = 1.0e-07;
  // U must be the global z axis 
  if( fabs(1.0 - V.dot(Z)) > eps ) {
    std::cout << "ILDPlanarMeasLayer: TrackerHitPlane measurment vectors V is not equal to the global Z axis. \n\n exit(1) called from file " << __FILE__ << " and line " << __LINE__ << std::endl;
    exit(1);
  }
  
  if( fabs(U.dot(Z)) > eps ) {
    std::cout << "ILDPlanarMeasLayer: TrackerHitPlane measurment vectors U is not in the global X-Y plane. \n\n exit(1) called from file " << __FILE__ << " and line " << __LINE__ << std::endl;
    exit(1);
  }

  // remember here the "position" of the hit in fact defines the origin of the plane it defines so u and v are per definition 0. 
  const edm4hep::Vector3d& pos=trkhit.getPosition();
  const TVector3 hit(pos.x, pos.y, pos.z) ;
  
  // convert to layer coordinates       
  TKalMatrix h    = this->XvToMv(hit);
  
  Double_t  x[2] ;
  Double_t dx[2] ;
  
  x[0] = h(0, 0);
  x[1] = h(1, 0);
  
  dx[0] = trkhit.getCovMatrix(2);
  dx[1] = trkhit.getCovMatrix(5);
  
  bool hit_on_surface = IsOnSurface(hit);
  /*
  std::cout << "ILDPlanarMeasLayer::ConvertLCIOTrkHit ILDPlanarHit created" 
			<< " for CellID " << trkhit.getCellID()
			<< " Layer R = " << this->GetXc().Perp() 
			<< " Layer phi = " << this->GetXc().Phi() 
			<< " Layer z0 = " << this->GetXc().Z() 
			<< " u = "  <<  x[0]
			<< " v = "  <<  x[1]
			<< " du = " << dx[0]
			<< " dv = " << dx[1]
			<< " x = " << hit.x()
			<< " y = " << hit.y()
			<< " z = " << hit.z()
			<< " r = " << hit.Perp()
			<< " onSurface = " << hit_on_surface
			<< std::endl ;
  */
  
  return hit_on_surface ? new ILDPlanarHit( *this , x, dx, this->GetBz(), trkhit) : NULL; 
}
