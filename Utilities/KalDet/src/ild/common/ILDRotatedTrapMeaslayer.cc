#include <iostream>

#include "ILDRotatedTrapMeaslayer.h"
#include "ILDPlanarHit.h"

#include "kaltest/TVTrack.h"
#include "TVector3.h"
#include "TMath.h"
#include "TRotMatrix.h"
#include "TBRIK.h"
#include "TNode.h"
#include "TString.h"

//#include <EVENT/TrackerHitPlane.h>

// #include "streamlog/streamlog.h"

ILDRotatedTrapMeaslayer::ILDRotatedTrapMeaslayer(TMaterial &min,
                                                 TMaterial &mout,
                                                 const TVector3  &center,
                                                 const TVector3  &normal,
                                                 Double_t   Bz,
                                                 Double_t   SortingPolicy,
                                                 Double_t   height,
                                                 Double_t   innerBaseLength,
                                                 Double_t   outerBaseLength,
                                                 Double_t   alpha,
                                                 Int_t      half_petal,
                                                 Bool_t     is_active,
                                                 Int_t      CellID,
                                                 const Char_t    *name)
: ILDVMeasLayer(min, mout, Bz, is_active, CellID, name),
TPlane(center, normal),
_sortingPolicy(SortingPolicy), _innerBaseLength(innerBaseLength), _outerBaseLength(outerBaseLength), _halfPetal(half_petal)
{
  
  if( GetXc().Z() >= 0 ) {
    _signZ = +1.0 ;
  }
  else {
    _signZ = -1.0 ;
  }
  
  _innerR = GetXc().Perp() - 0.5 * height ;
  _outerR = GetXc().Perp() + 0.5 * height ;
  _cosPhi = GetXc().X() / GetXc().Perp() ;
  _sinPhi = GetXc().Y() / GetXc().Perp() ;
  
  // alpha should be limited to +/- Pi/2
  _cosAlpha = cos(alpha) ;
  _sinAlpha = sin(alpha) ;
  
  _tanBeta = 0.5 * _outerBaseLength / _outerR  ;
  
}

TKalMatrix ILDRotatedTrapMeaslayer::XvToMv(const TVector3 &xv) const
{
  
  // Calculate measurement vector (hit coordinates) from global coordinates:
  
  TKalMatrix mv(kMdim,1);
  
  // SJA:FIXME: what to do with the -z hits, are they reflective, i.e. that means that the transverse coordinate with be reversed 
  // transverse coordinate, measured from the centre line of the petal ( the axis of symmetry of the trapizoid ) 
  mv(0,0)  = _cosAlpha * _signZ * ( _cosPhi*xv.Y() - _sinPhi*xv.X() ) + _sinAlpha * ( xv.Z() - GetXc().Z() ) ;
  
  // radial coordinate, measured from R = 0 ( x=0, y=0) 
  mv(1,0)  = _cosPhi * xv.X() + _sinPhi * xv.Y() ;
  return mv;
  
}


TVector3 ILDRotatedTrapMeaslayer::HitToXv(const TVTrackHit &vht) const
{
  //  const ILDPlanarHit &mv = dynamic_cast<const ILDPlanarHit &>(vht);
  
  Double_t x =   vht(1,0) * _cosPhi  - _signZ * _cosAlpha * vht(0,0) * _sinPhi  ;
  Double_t y =   vht(1,0) * _sinPhi  + _signZ * _cosAlpha * vht(0,0) * _cosPhi  ;
  
  Double_t z = GetXc().Z() + _signZ * vht(0,0) * _sinAlpha;
  
  return TVector3(x,y,z);
}

void ILDRotatedTrapMeaslayer::CalcDhDa(const TVTrackHit &vht,
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
  
  Int_t sdim = H.GetNcols();
  Int_t hdim = TMath::Max(5,sdim-1);
  
  // Set H = (@h/@a) = (@d/@a, @z/@a)^t
  
  for (Int_t i=0; i<hdim; i++) {
    
    H(0,i) = _cosAlpha * _signZ * ( _cosPhi*dxphiada(1,i) - _sinPhi*dxphiada(0,i) ) + _sinAlpha*dxphiada(2,i) ;
    H(1,i) = _cosPhi * dxphiada(0,i) + _sinPhi*dxphiada(1,i);
    
  }
  if (sdim == 6) {
    H(0,sdim-1) = 0.;
    H(1,sdim-1) = 0.;
  }
  
}



Bool_t ILDRotatedTrapMeaslayer::IsOnSurface(const TVector3 &xx) const
{
  
  //  streamlog_out(DEBUG0) << "IsOnSurface " << std::endl;  
  
  bool onSurface = false ;
  
  TKalMatrix mv = XvToMv(xx);
  
  // check whether the hit lies in the same plane as the surface
  if( TMath::Abs((xx.X()-GetXc().X())*GetNormal().X() + (xx.Y()-GetXc().Y())*GetNormal().Y() + (xx.Z()-GetXc().Z())*GetNormal().Z()) < 1e-4){
    // check whether the hit lies within the boundary of the surface 
    if(  mv(1,0) <= _outerR   &&  mv(1,0) >= _innerR 
       && 
       TMath::Abs(mv(0,0)) <=   mv(1,0) * _tanBeta  )
        { 
          
          // meaning of _halfPetal:
          //                  0 complete trapezoid
          //                 +1 positive half only, i.e. the side of the petal in which the transverse coordinate, mv(0,0), is positive 
          //                 -1 negative half only, i.e. the side of the petal in which the transverse coordinate, mv(0,0), is negative
          
          if( _halfPetal == 0 || ( _halfPetal * mv(0,0) ) >= 0) { // check if the point lies in the correct half
            onSurface = true ;
          }
          
        }
    
    else{
      onSurface = false;
    }
    
  }
  
  return onSurface;
  
}


ILDVTrackHit* ILDRotatedTrapMeaslayer::ConvertLCIOTrkHit(edm4hep::ConstTrackerHit trkhit) const {
  
  //EVENT::TrackerHitPlane* plane_hit = dynamic_cast<EVENT::TrackerHitPlane*>( trkhit ) ;
  if((trkhit.getType()&8)!=8) return NULL;
  //if( plane_hit == NULL )  return NULL; // SJA:FIXME: should be replaced with an exception  
  const edm4hep::Vector3d& pos=trkhit.getPosition();
  const TVector3 hit(pos.x, pos.y, pos.z);
  
  // convert to layer coordinates       
  TKalMatrix h    = this->XvToMv(hit);
  
  Double_t  x[2] ;
  Double_t dx[2] ;
  
  x[0] = h(0, 0);
  x[1] = h(1, 0);
  
  dx[0] = trkhit.getCovMatrix(2);
  dx[1] = trkhit.getCovMatrix(5);
  
  bool hit_on_surface = IsOnSurface(hit);
  
  // streamlog_out(DEBUG1) << "ILDRotatedTrapMeaslayer::ConvertLCIOTrkHit ILDPlanarHit created" 
  //       		<< " u = "  <<  x[0]
  //       		<< " v = "  <<  x[1]
  //       		<< " du = " << dx[0]
  //       		<< " dv = " << dx[1]
  //       		<< " x = " << hit.x()
  //       		<< " y = " << hit.y()
  //       		<< " z = " << hit.z()
  //       		<< " onSurface = " << hit_on_surface
  //       		<< std::endl ;
  
  return hit_on_surface ? new ILDPlanarHit( *this , x, dx, this->GetBz(),trkhit) : NULL; 
  
  
}
