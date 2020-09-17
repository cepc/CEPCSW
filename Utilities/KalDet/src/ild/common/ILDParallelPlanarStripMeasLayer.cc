
#include "ILDParallelPlanarStripMeasLayer.h"
#include "ILDPlanarStripHit.h"
#include "kaltest/TVTrack.h"
#include "kaltest/TVTrackHit.h"

//#include <EVENT/TrackerHitPlane.h>

#include "gearimpl/Vector3D.h"

// #include "streamlog/streamlog.h"


TKalMatrix ILDParallelPlanarStripMeasLayer::XvToMv(const TVector3 &xv) const
{
 
//  streamlog_out(DEBUG0) << "\t ILDParallelPlanarStripMeasLayer::XvToMv: "
//  << " x = " << xv.X() 
//  << " y = " << xv.Y() 
//  << " z = " << xv.Z() 
//  << std::endl;

  
  double cos_phi = GetNormal().X()/GetNormal().Perp();
  double sin_phi = GetNormal().Y()/GetNormal().Perp();
  
  // theta is the strip angle rotation in the plane of the wafer
  double cos_theta = cos(_stripAngle);
  double sin_theta = sin(_stripAngle);
  
  double delta_x = xv.X() - GetXc().X();
  double delta_y = xv.Y() - GetXc().Y();

  double delta_t = (delta_x * sin_phi - delta_y * cos_phi) ; 
  double delta_z = xv.Z() - GetXc().Z();


  TKalMatrix mv(ILDPlanarStripHit_DIM,1);
  
  mv(0,0) = ( (delta_t + fUOrigin) * cos_theta ) + delta_z * sin_theta ;
  
  if (ILDPlanarStripHit_DIM == 2) {
    mv(1,0) = delta_z * cos_theta - delta_t * sin_theta;
  }
  
//  streamlog_out(DEBUG0) << "\t ILDParallelPlanarStripMeasLayer::XvToMv: " 
//  << " mv(0,0) = " << mv(0,0) ;
//  if (ILDPlanarStripHit_DIM == 2) {
//    streamlog_out(DEBUG0) << " mv(1,0) = " << mv(1,0);
//  }
//  streamlog_out(DEBUG0) << std::endl;

  
  return mv;

}


TVector3 ILDParallelPlanarStripMeasLayer::HitToXv(const TVTrackHit &vht) const
{

//  streamlog_out(DEBUG0) << "\t ILDParallelPlanarStripMeasLayer::HitToXv: "
//  << " vht(0,0) = " << vht(0,0);
//  if (ILDPlanarStripHit_DIM == 2) {
//    streamlog_out(DEBUG0) << " vht(1,0) = " << vht(1,0);
//  }
//  streamlog_out(DEBUG0) << std::endl;
  
  double cos_phi = GetNormal().X()/GetNormal().Perp();
  double sin_phi = GetNormal().Y()/GetNormal().Perp();
  
  // theta is the strip angle rotation in the plane of the wafer
  double cos_theta = cos(_stripAngle);
  double sin_theta = sin(_stripAngle);
  
  double t =  (vht(0,0) - fUOrigin ) * cos_theta ;
    
  double x =  t * sin_phi + this->GetXc().X();
  double y = -t * cos_phi + this->GetXc().Y();
  
  double dz = 0.0;
  
  if (ILDPlanarStripHit_DIM == 2) {
    dz = vht(1,0) * cos_theta ;
  }

  
  double z = vht(0,0) * sin_theta + this->GetXc().Z() + dz;

//  streamlog_out(DEBUG0) << "\t ILDParallelPlanarStripMeasLayer::HitToXv: "
//  << " x = " << x 
//  << " y = " << y 
//  << " z = " << z 
//  << std::endl;
  
  
  this->IsOnSurface(TVector3(x,y,z));
  
  return TVector3(x,y,z);

}

void ILDParallelPlanarStripMeasLayer::CalcDhDa(const TVTrackHit &vht,
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

  // theta is the strip angle rotation in the plane of the wafer
  double cos_theta = cos(_stripAngle);
  double sin_theta = sin(_stripAngle);
  
  Int_t sdim = H.GetNcols();
  Int_t hdim = TMath::Max(5,sdim-1);
  
  // Set H = (@h/@a) = (@d/@a, @z/@a)^t

  double dudx =  cos_theta * sin_phi;
  double dudy = -cos_theta * cos_phi;
  double dudz =  sin_theta;

  double dvdx = -sin_theta * sin_phi;
  double dvdy =  sin_theta * cos_phi;
  double dvdz =  cos_theta;

  
//  streamlog_out(DEBUG0) << "\t ILDParallelPlanarStripMeasLayer::CalcDhDa: dudx = " << dudx << " dudy = " << dudy << " dudz = " << dudz << std::endl ;
//  if (ILDPlanarStripHit_DIM == 2) streamlog_out(DEBUG0) << "\t ILDParallelPlanarStripMeasLayer::CalcDhDa: dvdx = " << dvdx << " dvdy = " << dvdy << " dvdz = " << dvdz << std::endl ;
  
  for (Int_t i=0; i<hdim; i++) {
    
    H(0,i) =  dudx * dxphiada(0,i) + dudy * dxphiada(1,i) + dudz * dxphiada(2,i) ;   

//    streamlog_out(DEBUG0)  << " H(0,"<< i <<") = " << H(0,i)  ;

    if (ILDPlanarStripHit_DIM == 2) {

      H(1,i) = dvdx * dxphiada(0,i) + dvdy * dxphiada(1,i) + dvdz * dxphiada(2,i);
//      streamlog_out(DEBUG0)  << " H(1,"<< i <<") = " << H(1,i);

    }
    
  }
  if (sdim == 6) {
    H(0,sdim-1) = 0.;
    if (ILDPlanarStripHit_DIM == 2) {
      H(1,sdim-1) = 0.;
    }
  }  
  
//  streamlog_out(DEBUG0) << std::endl;
  
}

ILDVTrackHit* ILDParallelPlanarStripMeasLayer::ConvertLCIOTrkHit(edm4hep::ConstTrackerHit trkhit) const {
  
  //EVENT::TrackerHitPlane* plane_hit = dynamic_cast<EVENT::TrackerHitPlane*>( trkhit ) ;
  if((trkhit.getType()&8)!=8) {
    //if( plane_hit == NULL )  { 
    std::cout << "ILDParallelPlanarStripMeasLayer::ConvertLCIOTrkHit dynamic_cast to ILDPlanarStripHit failed " << std::endl; 
    throw std::logic_error("Invalid invoke ILDParallelPlanarStripMeasLayer by TrackerHit trkhit");
  }
  
  
  // remember here the "position" of the hit in fact defines the origin of the plane it defines so u and v are per definition 0. 
  // this is still the case for a 1-dimentional measurement, and is then used to calculate the u coordinate according to the origin of the actual measurement plane.
  //const TVector3 hit( plane_hit.getPosition()[0], plane_hit.getPosition()[1], plane_hit.getPosition()[2]) ;
  const edm4hep::Vector3d& pos=trkhit.getPosition(); 
  const TVector3 hit(pos.x, pos.y, pos.z);
  
  // convert to layer coordinates       
  TKalMatrix h(ILDPlanarStripHit_DIM,1);    
  h = this->XvToMv(hit);
  
  Double_t  x[ILDPlanarStripHit_DIM] ;
  Double_t dx[ILDPlanarStripHit_DIM] ;
  
  x[0] = h(0, 0);
  if(ILDPlanarStripHit_DIM == 2) x[1] = h(1, 0);

  //dx[0] = plane_hit.getdU() ;
  //if(ILDPlanarStripHit_DIM == 2) dx[1] = plane_hit.getdV() ;
  dx[0] = trkhit.getCovMatrix(2);
  if(ILDPlanarStripHit_DIM == 2) dx[1] = trkhit.getCovMatrix(5);

  bool hit_on_surface = IsOnSurface(hit);
  /*
  std::cout << "ILDParallelPlanarStripMeasLayer::ConvertLCIOTrkHit ILDPlanarStripHit created" 
	    << " for CellID " << trkhit.getCellID()
	    << " Layer R = " << this->GetXc().Perp() 
	    << " Layer phi = " << this->GetXc().Phi() 
	    << " Layer z0 = " << this->GetXc().Z() 
	    << " u = "  <<  x[0]
	    << " du = " << dx[0];
	    
  if(ILDPlanarStripHit_DIM == 2)  std::cout << " v = "  <<  x[1] << " dv = " << dx[1];
  
  std::cout << " x = " << pos.x
	    << " y = " << pos.y
	    << " z = " << pos.z
	    << " r = " << sqrt(pos.x*pos.x + pos.y*pos.y)
	    << " onSurface = " << hit_on_surface
	    << std::endl ;
  */
  ILDPlanarStripHit hh( *this , x, dx, this->GetBz(),trkhit);
  
  this->HitToXv(hh);
  
  return hit_on_surface ? new ILDPlanarStripHit( *this , x, dx, this->GetBz(),trkhit) : NULL; 
  
  
}


