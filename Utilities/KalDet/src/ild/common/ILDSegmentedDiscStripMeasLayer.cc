
#include "ILDSegmentedDiscStripMeasLayer.h"

#include "ILDPlanarStripHit.h"

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

#include <math.h>
#include <assert.h>
#include <algorithm>

// #include "streamlog/streamlog.h"


TKalMatrix ILDSegmentedDiscStripMeasLayer::XvToMv(const TVector3 &xv) const
{
  
  // Calculate measurement vector (hit coordinates) from global coordinates:
  
  // streamlog_out(DEBUG0) << "\t ILDSegmentedDiscStripMeasLayer::XvToMv: "
  // << " x = " << xv.X() 
  // << " y = " << xv.Y() 
  // << " z = " << xv.Z() 
  // << std::endl;
  
  // let's start with the sensor whose axis of symmetry is 
  // aligned with the y-axis and whose sensitive face is facing towards the IP.
  // For a zero strip angle then:
  //                              u = -x
  //                              v =  y
  //                              w = -z
  
  // coordinate matrix to return
  TKalMatrix mv(ILDPlanarStripHit_DIM,1);
  
  int segmentIndex = get_segment_index(xv.Phi()) ;
  
  TVector3 XC = this->get_segment_centre(segmentIndex);
  
  double sensor_x0 = XC.X();
  double sensor_y0 = XC.Y();
  
  // here we are assuming that there is no offset of the centre of the sensor in the x-y plane. 
  // SJA:FIXME: We need to get the segment we are in to get phi
  
  double phi_sensor = XC.Phi();
  
  // for ILDSegmentedDiscMeasLayer the Normal is pointing from the IP to the Plane
  double sign_z = GetNormal().Z() < 0 ? -1.0 : 1.0 ;

  double phi = phi_sensor + sign_z*M_PI/2 ;
  
  phi += _stripAngle; // the strip rotation is around the w vector which is pointing at the IP
  
  double cos_phi = cos(phi);
  double sin_phi = sin(phi);

  double delta_x = xv.X() - sensor_x0;
  double delta_y = xv.Y() - sensor_y0;
  
  // is the  sensitive facing the IP (+1) or facing away from it (-1)
  //  const double sign_z = GetNormal().Z() < 0 ? 1.0 : -1.0 ;

  //  double u = sign_z * (delta_y * cos_phi - delta_x * sin_phi) ; 

  double cos_theta = -sign_z;
  
  //  double u = delta_x * cos_phi - cos_theta * delta_y * sin_phi ; 
  double u = delta_x * cos_phi + delta_y * sin_phi ; 

  mv(0,0) = u ;

  if (ILDPlanarStripHit_DIM == 2) {
    //    double v =   (delta_x * sin_phi + cos_theta * delta_y * cos_phi) ; 
    double v =   ( cos_theta * delta_y * cos_phi - cos_theta * delta_x * sin_phi) ; 
    mv(1,0)  = v ;
  }

  
  // streamlog_out(DEBUG0) << "\t ILDSegmentedDiscStripMeasLayer::XvToMv: phi_sensor = " << phi_sensor << " phi = " << phi << " stripAngle = " << _stripAngle << " sign_z = " << sign_z<< std::endl;
  
  // streamlog_out(DEBUG0) << "\t ILDSegmentedDiscStripMeasLayer::XvToMv: " 
  // << " mv(0,0) = " << mv(0,0) ;
  // if (ILDPlanarStripHit_DIM == 2) {
  //   streamlog_out(DEBUG0) << " mv(1,0) = " << mv(1,0);
  // }
  // streamlog_out(DEBUG0) << std::endl;
  
  return mv;
  
}

TVector3 ILDSegmentedDiscStripMeasLayer::HitToXv(const TVTrackHit &vht) const
{
  
  // streamlog_out(DEBUG0) << "\t ILDSegmentedDiscStripMeasLayer::HitToXv: "
  // << " vht(0,0) = " << vht(0,0);
  // if (ILDPlanarStripHit_DIM == 2) {
  //   streamlog_out(DEBUG0) << " vht(1,0) = " << vht(1,0);
  // }
  // streamlog_out(DEBUG0) << std::endl;

  
  const ILDPlanarStripHit &mv = dynamic_cast<const ILDPlanarStripHit &>(vht);
  
  UTIL::BitField64 encoder( lcio::ILDCellID0::encoder_string ) ;
  edm4hep::ConstTrackerHit hit = mv.getLCIOTrackerHit();
  encoder.setValue(hit.getCellID());
  int segmentIndex = encoder[lcio::ILDCellID0::module] / 2 ;
  

  TVector3 XC = this->get_segment_centre(segmentIndex);
  
  double sensor_x0 = XC.X();
  double sensor_y0 = XC.Y();
  double sensor_z0 = XC.Z();

//  streamlog_out(DEBUG0) << "\t ILDSegmentedDiscStripMeasLayer::HitToXv: segmentIndex = " << segmentIndex << " x0 = " << sensor_x0 << " y0 = " << sensor_y0 << " z0 = " << sensor_z0 << " segment Phi = " << XC.Phi() << std::endl; 
  
  // here we are assuming that there is no offset of the centre of the sensor in the x-y plane. 
  // SJA:FIXME: We need to get the segment we are in to get phi

  double phi_sensor = XC.Phi();
  // for ILDSegmentedDiscMeasLayer the Normal is pointing from the IP to the Plane
  double sign_z = GetNormal().Z() < 0 ? -1.0 : 1.0 ;;
  
  double phi = phi_sensor + sign_z*M_PI/2 ;
  
  phi += _stripAngle; // the strip rotation is around the w vector which is pointing at the IP
  
  double cos_phi = cos(phi);
  double sin_phi = sin(phi);

  double cos_theta = -sign_z;
    
  double u = mv(0,0);
  double v = 0.0;
  
  if (ILDPlanarStripHit_DIM == 2) {
    v = mv(1,0);
  }
  
//  double delta_x =  (u * cos_phi - cos_theta * v * sin_phi) ;
//  double delta_y =  (u * sin_phi + cos_theta * v * cos_phi) ;
  
  double delta_x =  (u * cos_phi - cos_theta * v * sin_phi) ;
  double delta_y =  (u * sin_phi + cos_theta * v * cos_phi) ;
  
  double x = delta_x + sensor_x0; 
  double y = delta_y + sensor_y0; 
  
  double z = sensor_z0 ;
  
  // streamlog_out(DEBUG0) << "\t ILDSegmentedDiscStripMeasLayer::HitToXv: "
  // << " x = " << x 
  // << " y = " << y 
  // << " z = " << z 
  // << std::endl;

  
  return TVector3(x,y,z);
}

void ILDSegmentedDiscStripMeasLayer::CalcDhDa(const TVTrackHit &vht,
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
  
  int segmentIndex = get_segment_index(xxv.Phi()) ;
  
  TVector3 XC = this->get_segment_centre(segmentIndex);

  
  double phi_sensor = XC.Phi();
  // for ILDSegmentedDiscMeasLayer the Normal is pointing from the IP to the Plane 
  double sign_z = GetNormal().Z() < 0 ? -1.0 : 1.0 ;
  
  double phi = phi_sensor + sign_z*M_PI/2 ;
  
  phi += _stripAngle; // the strip rotation is around the w vector which is pointing at the IP
  
  double cos_phi = cos(phi);
  double sin_phi = sin(phi);

  double cos_theta = -sign_z;
  
  Int_t sdim = H.GetNcols();
  Int_t hdim = TMath::Max(5,sdim-1);
  
  // Set H = (@h/@a) = (@d/@a, @z/@a)^t


  double dudx =   cos_phi;
  double dudy =   sin_phi;
  
  double dvdx =  -cos_theta * sin_phi;
  double dvdy =   cos_theta * cos_phi;

  // streamlog_out(DEBUG0) << "\t ILDSegmentedDiscStripMeasLayer::CalcDhDa: "
  // << " dudx = " << dudx 
  // << " dudy = " << dudy
  // << " dvdx = " << dvdx 
  // << " dvdy = " << dvdy 
  // << std::endl;
   
  for (Int_t i=0; i<hdim; i++) {
    
    H(0,i) = dudx * dxphiada(0,i) + dudy * dxphiada(1,i) ;
    
    if (ILDPlanarStripHit_DIM == 2) {
      H(1,i) = dvdx * dxphiada(0,i) + dvdy * dxphiada(1,i) ;
    }
    
  }
  if (sdim == 6) {
    H(0,sdim-1) = 0.0;
    if (ILDPlanarStripHit_DIM == 2) {
      H(1,sdim-1) = 0.;
    }

  }
  
}



ILDVTrackHit* ILDSegmentedDiscStripMeasLayer::ConvertLCIOTrkHit(edm4hep::ConstTrackerHit trkhit) const {
  
  //EVENT::TrackerHitPlane* plane_hit = dynamic_cast<EVENT::TrackerHitPlane*>( trkhit ) ;
  if((trkhit.getType()&8)!=8){
  //if( plane_hit == NULL )  { 
    // streamlog_out(ERROR) << "ILDSegmentedDiscStripMeasLayer::ConvertLCIOTrkHit dynamic_cast to TrackerHitPlane failed " << std::endl; 
    return NULL; // SJA:FIXME: should be replaced with an exception  
  }
  
  // remember here the "position" of the hit in fact defines the origin of the plane it defines so u and v are per definition 0. 
  // this is still the case for a 1-dimentional measurement, and is then used to calculate the u coordinate according to the origin of the actual measurement plane.
  const edm4hep::Vector3d& pos=trkhit.getPosition();
  const TVector3 hit(pos[0], pos[1], pos[2]) ;
  
  // convert to layer coordinates       
  TKalMatrix h(ILDPlanarStripHit_DIM,1);  
  
  h = this->XvToMv(hit);
  
  double  x[ILDPlanarStripHit_DIM] ;
  double dx[ILDPlanarStripHit_DIM] ;
  
  x[0] = h(0, 0);
  if(ILDPlanarStripHit_DIM == 2) x[1] = h(1, 0);
  
  dx[0] = trkhit.getCovMatrix(2);
  if(ILDPlanarStripHit_DIM == 2) dx[1] = trkhit.getCovMatrix(5);
    
  bool hit_on_surface = IsOnSurface(hit);
  
  //std::cout << "ILDSegmentedDiscStripMeasLayer::ConvertLCIOTrkHit ILDPlanarStripHit created" 
  //	    << " for CellID " << trkhit.getCellID()
  //	    << " Disc Z = " << this->GetXc().Z() 
  //	    << " u = "  <<  x[0]
  //	    << " du = " << dx[0];
  
  //if(ILDPlanarStripHit_DIM == 2)  std::cout << " v = "  <<  x[1] << " dv = " << dx[1];
  
  //std::cout << " x = " << hit.x()
  //	    << " y = " << hit.y()
  //	    << " z = " << hit.z()
  //	    << " r = " << hit.Perp()
  //	    << " onSurface = " << hit_on_surface
  //	    << std::endl ;
  
  ILDPlanarStripHit hh( *this , x, dx, this->GetBz(),trkhit);
  
  this->HitToXv(hh);
  
  return hit_on_surface ? new ILDPlanarStripHit( *this , x, dx, this->GetBz(),trkhit) : NULL; 
  
}



