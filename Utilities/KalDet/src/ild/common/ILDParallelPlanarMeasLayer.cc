
#include "ILDParallelPlanarMeasLayer.h"
#include "kaltest/TVTrack.h"

// #include "streamlog/streamlog.h"


Int_t ILDParallelPlanarMeasLayer::CalcXingPointWith(const TVTrack  &hel,
                                                    TVector3 &xx,
                                                    Double_t &phi,
                                                    Int_t     mode,
                                                    Double_t  eps) const{
  
  
  // check that direction has one of the correct values
  if( !( mode == 0 || mode == 1 || mode == -1) ) return -1 ;
  
  
  
  // This assumes nonzero B field.
  //
  // Copy helix parameters to local variables.
  //
  
  Double_t dr     = hel.GetDrho();
  Double_t phi0   = hel.GetPhi0(); //
  Double_t kappa  = hel.GetKappa();
  Double_t rho    = hel.GetRho();
  Double_t omega  = 1.0 / rho;
  Double_t r      = TMath::Abs(rho);
  Double_t z0     = hel.GetDz();
  Double_t tanl   = hel.GetTanLambda();
  
  TVector3 ref_point = hel.GetPivot();
  
  //
  // Check if charge is nonzero.
  //
  
  Int_t    chg = (Int_t)TMath::Sign(1.1,kappa);
  if (!chg) {
    // streamlog_out(ERROR) << ">>>> Error >>>> ILDParallelPlanarMeasLayer::CalcXingPointWith" << std::endl
    // << "      Kappa = 0 is invalid for a helix "          << std::endl;
    return -1;
  }
  
  //
  // Project everything to XY plane and calculate crossing points.
  //
  // taken from http://paulbourke.net/geometry/sphereline/
  
  
  const double sin_phi0 = sin(phi0); 
  const double cos_phi0 = cos(phi0); 
  
  const double x_pca = ref_point.x() + dr * cos_phi0 ; 
  const double y_pca = ref_point.y() + dr * sin_phi0 ; 
  const double z_pca = ref_point.z() + z0 ;
  
  const double x_c   = ref_point.x() + ( rho + dr) * cos_phi0 ;
  const double y_c   = ref_point.y() + ( rho + dr) * sin_phi0 ;
  
  // get the extreams of the plane in x and y
  const double x1 = this->GetXc().x() - 0.5*this->GetXiwidth()*sin(this->GetXc().Phi()) - this->GetXioffset() * sin(this->GetXc().Phi());
  const double x2 = this->GetXc().x() + 0.5*this->GetXiwidth()*sin(this->GetXc().Phi()) - this->GetXioffset() * sin(this->GetXc().Phi());
  
  const double y1 = this->GetXc().y() + 0.5*this->GetXiwidth()*cos(this->GetXc().Phi()) + this->GetXioffset() * cos(this->GetXc().Phi());
  const double y2 = this->GetXc().y() - 0.5*this->GetXiwidth()*cos(this->GetXc().Phi()) + this->GetXioffset() * cos(this->GetXc().Phi());
  
  //  streamlog_out(DEBUG0) << "ILDParallelPlanarMeasLayer::CalcXingPointWith: Xc = " << this->GetXc().x() << " Yc = " << this->GetXc().y() << " Rc = " << this->GetXc().Perp() << " Phi = " << this->GetXc().Phi() << std::endl;
  //    streamlog_out(DEBUG0) << "ILDParallelPlanarMeasLayer::CalcXingPointWith: x1 = " << x1 << " y1 = "  << y1 << " R " << TVector3(x1,y1,0).Perp() << " phi = " << TVector3(x1,y1,0).Phi() << std::endl ; 
  //    streamlog_out(DEBUG0) << "ILDParallelPlanarMeasLayer::CalcXingPointWith: x2 = " << x2 << " y2 = "  << y2 << " R " << TVector3(x2,y2,0).Perp() << " phi = " << TVector3(x2,y2,0).Phi() << std::endl ; 
  
  const double dx = x2 - x1 ;
  const double dy = y2 - y1 ;
  
  const double a = dx * dx + dy * dy ;
  
  const double b = 2.0 * ( dx * (x1 - x_c) + dy * (y1 - y_c) ) ;
  
  double c = x_c * x_c + y_c * y_c;
  c += x1 * x1 + y1 * y1 ;
  
  c -= 2.0 * ( x_c * x1 + y_c * y1 );
  c -= rho * rho;
  
  const double bb4ac = b * b - 4.0 * a * c;
  
  double u1 ; // first solution
  double u2 ; // second solution
  
  double x_ins = DBL_MAX; 
  double y_ins = DBL_MAX;
  double s_ins = DBL_MAX;
  
  /* Check for solvability. */
  if (bb4ac + eps < 0.0 ) { // no intersection
    return 0;
  }
  else if(bb4ac - eps < 0.0) { // circle intersects at one point, tangential 
    return 0;
  }
  else{
    u1 = (-b + sqrt(bb4ac)) / (2.0 * a);
    u2 = (-b - sqrt(bb4ac)) / (2.0 * a);
  }
  
  // test values of u 
  // case i)   ( u1 < 0 && u2 < 0 ) || ( u1 > 1 && u2 > 1 )
  // Line segment doesn't intersect and is on outside of the circle 
  if ( ( u1 <= 0 && u2 <= 0 ) || ( u1 >= 1 && u2 >= 1 ) ) {
    
    //          streamlog_out(DEBUG0) << "ILDParallelPlanarMeasLayer::CalcXingPointWith: Line segment doesn't intersect and is outside the circle " << std::endl;
    //          const double x_ins1 = x1 + u1 * (dx) ;
    //          const double y_ins1 = y1 + u1 * (dy) ;
    //          
    //          const double x_ins2 = x1 + u2 * (dx) ;
    //          const double y_ins2 = y1 + u2 * (dy) ;
    //          
    //          streamlog_out(DEBUG0) << "ILDParallelPlanarMeasLayer::CalcXingPointWith: 1st solution u = " << u1 << " : " << u1 * dx << " " << u1 * dy << std::endl ;  
    //          streamlog_out(DEBUG0) << "ILDParallelPlanarMeasLayer::CalcXingPointWith: 2nd solution u = " << u2 << " : " << u2 * dx << " " << u2 * dy << std::endl ;
    //          streamlog_out(DEBUG0) << "ILDParallelPlanarMeasLayer::CalcXingPointWith: 1st solution x:y = " << x_ins1 << " : " << y_ins1 << std::endl ;
    //          streamlog_out(DEBUG0) << "ILDParallelPlanarMeasLayer::CalcXingPointWith: 2nd solution x:y = " << x_ins2 << " : " << y_ins2 << std::endl ;
    
    return 0 ;
  }
  
  // case ii)  ( u1 < 0 && u2 > 1 ) || ( u1 > 1 && u2 < 0 )
  // Line segment doesn't intersect and is inside the circle 
  else if ( ( u1 <= 0 && u2 >= 1 ) || ( u1 >= 1 && u2 <= 0 ) ) {
    
    //    streamlog_out(DEBUG0) << "ILDParallelPlanarMeasLayer::CalcXingPointWith: Line segment doesn't intersect and is inside the circle " << std::endl;
    //          
    //          const double x_ins1 = x1 + u1 * (dx) ;
    //          const double y_ins1 = y1 + u1 * (dy) ;
    //          
    //          const double x_ins2 = x1 + u2 * (dx) ;
    //          const double y_ins2 = y1 + u2 * (dy) ;
    //          
    //          streamlog_out(DEBUG0) << "ILDParallelPlanarMeasLayer::CalcXingPointWith: 1st solution u = " << u1 << " : " << u1 * dx << " " << u1 * dy << std::endl ;  
    //          streamlog_out(DEBUG0) << "ILDParallelPlanarMeasLayer::CalcXingPointWith: 2nd solution u = " << u2 << " : " << u2 * dx << " " << u2 * dy << std::endl ;
    //          streamlog_out(DEBUG0) << "ILDParallelPlanarMeasLayer::CalcXingPointWith: 1st solution x:y = " << x_ins1 << " : " << y_ins1 << std::endl ;
    //          streamlog_out(DEBUG0) << "ILDParallelPlanarMeasLayer::CalcXingPointWith: 2nd solution x:y = " << x_ins2 << " : " << y_ins2 << std::endl ;
    
    
    return 0 ;
  }
  
  // case iii) ( u1 > 0 && u1 < 1 && u2 < 0 && u2 > 1 ) || ( u2 > 0 && u2 < 1 && u1 < 0 && u1 > 1 )
  // Line segment intersects at one point
  else if ( ( u1 > 0 && u1 < 1 && (u2 <= 0 || u2 >= 1) ) || ( u2 > 0 && u2 < 1 && ( u1 <= 0 || u1 >= 1 )) ) {
    
    //          streamlog_out(DEBUG0) << "ILDParallelPlanarMeasLayer::CalcXingPointWith: Only one possible solution" << std::endl;
    //          streamlog_out(DEBUG0) << "ILDParallelPlanarMeasLayer::CalcXingPointWith: 1st solution u = " << u1 << " : " << u1 * dx << " " << u1 * dy << std::endl ;  
    //          streamlog_out(DEBUG0) << "ILDParallelPlanarMeasLayer::CalcXingPointWith: 2nd solution u = " << u2 << " : " << u2 * dx << " " << u2 * dy << std::endl ;
    
    if ( u1 > 0 && u1 < 1 ) { // use u1 
      
      x_ins = x1 + u1 * (dx) ;
      y_ins = y1 + u1 * (dy) ;
      
      //                        streamlog_out(DEBUG0) << "ILDParallelPlanarMeasLayer::CalcXingPointWith: take 1st solution x:y:z = " << x_ins << " : " << y_ins ;
    }
    else{ // use u2
      
      x_ins = x1 + u2 * (dx) ;
      y_ins = y1 + u2 * (dy) ;
      
      //                        streamlog_out(DEBUG0) << "ILDParallelPlanarMeasLayer::CalcXingPointWith: take 2nd solution x:y:z = " << x_ins << " : " << y_ins ;
    }
    
    const double delta_x = x_ins - x_pca ;
    const double delta_y = y_ins - y_pca ;
    
    const double sin_delta_phi =       omega*delta_x*sin_phi0 - omega*delta_y*cos_phi0 ;
    const double cos_delta_phi = 1.0 - omega*delta_x*cos_phi0 - omega*delta_y*sin_phi0 ;
    
    //          const double sin_delta_phi =     - omega*delta_x*cos_phi0 - omega*delta_y*sin_phi0 ;
    //    const double cos_delta_phi = 1.0 - omega*delta_x*sin_phi0 + omega*delta_y*cos_phi0 ;
    
    s_ins = atan2(-sin_delta_phi,cos_delta_phi) / omega ;
    
    //SJA:FIXME: do we need to consider the mode here ...
    
    //          streamlog_out(DEBUG0) << " : " << z_pca + s_ins * tanl << std::endl ;
    
  }
  
  // case iv)  ( u1 > 0 && u1 < 1 && u2 > 0 && u2 < 1 ) 
  // Line segment intersects at two points
  else{
    
    const double x_ins1 = x1 + u1 * (dx) ;
    const double y_ins1 = y1 + u1 * (dy) ;
    
    const double x_ins2 = x1 + u2 * (dx) ;
    const double y_ins2 = y1 + u2 * (dy) ;
    
    //          streamlog_out(DEBUG0) << "ILDParallelPlanarMeasLayer::CalcXingPointWith: 1st solution u = " << u1 << " : " << u1 * dx << " " << u1 * dy << std::endl ;  
    //          streamlog_out(DEBUG0) << "ILDParallelPlanarMeasLayer::CalcXingPointWith: 2nd solution u = " << u2 << " : " << u2 * dx << " " << u2 * dy << std::endl ;
    //          streamlog_out(DEBUG0) << "ILDParallelPlanarMeasLayer::CalcXingPointWith: 1st solution x:y = " << x_ins1 << " : " << y_ins1 << std::endl ;
    //          streamlog_out(DEBUG0) << "ILDParallelPlanarMeasLayer::CalcXingPointWith: 2nd solution x:y = " << x_ins2 << " : " << y_ins2 << std::endl ;
    
    // now calculate the path lengths 
    double s_1 = 0.0 ;
    double s_2 = 0.0 ;
    
    const double delta_x1 = x_ins1 - x_pca ;
    const double delta_y1 = y_ins1 - y_pca ;
    
    //          const double sin_delta_phi1 =     - omega*delta_x1*cos_phi0 - omega*delta_y1*sin_phi0 ;
    //    const double cos_delta_phi1 = 1.0 - omega*delta_x1*sin_phi0 + omega*delta_y1*cos_phi0 ;
    
    const double sin_delta_phi1 =       omega*delta_x1*sin_phi0 - omega*delta_y1*cos_phi0 ;
    const double cos_delta_phi1 = 1.0 - omega*delta_x1*cos_phi0 - omega*delta_y1*sin_phi0 ;
    
    s_1 = atan2(-sin_delta_phi1,cos_delta_phi1) / omega ;
    
    
    const double delta_x2 = x_ins2 - x_pca ;
    const double delta_y2 = y_ins2 - y_pca ;
    
    //          const double sin_delta_phi2 =     - omega*delta_x2*cos_phi0 - omega*delta_y2*sin_phi0 ;
    //    const double cos_delta_phi2 = 1.0 - omega*delta_x2*sin_phi0 + omega*delta_y2*cos_phi0 ;
    
    const double sin_delta_phi2 =       omega*delta_x2*sin_phi0 - omega*delta_y2*cos_phi0 ;
    const double cos_delta_phi2 = 1.0 - omega*delta_x2*cos_phi0 - omega*delta_y2*sin_phi0 ;
    
    s_2 = atan2(-sin_delta_phi2,cos_delta_phi2) / omega ;
    
    
    if( mode == 0 ) { // take closest intersection
      if( TMath::Abs(s_1) < TMath::Abs(s_2) ) {
        //                              streamlog_out(DEBUG0) << "ILDParallelPlanarMeasLayer::CalcXingPointWith: take 1st solution " << std::endl;
        x_ins = x_ins1;
        y_ins = y_ins1;
        s_ins = s_1;
      }
      else {
        //                              streamlog_out(DEBUG0) << "ILDParallelPlanarMeasLayer::CalcXingPointWith: take 2nd solution " << std::endl;
        x_ins = x_ins2;
        y_ins = y_ins2;
        s_ins = s_2;
      }
    }
    
    else{
      
      if ( s_1 < 0.0 ) s_1 +=  TMath::TwoPi() * r ;
      if ( s_2 < 0.0 ) s_2 +=  TMath::TwoPi() * r ;
      
      if( mode == 1 ){ // take the intersection with smallest s
        if( s_1 < s_2 ) {
          //                                    streamlog_out(DEBUG0) << "ILDParallelPlanarMeasLayer::CalcXingPointWith: take 1st solution " << std::endl;
          x_ins = x_ins1;
          y_ins = y_ins1;
          s_ins = s_1;
        }
        else {
          //                                    streamlog_out(DEBUG0) << "ILDParallelPlanarMeasLayer::CalcXingPointWith: take 2nd solution " << std::endl;
          x_ins = x_ins2;
          y_ins = y_ins2;
          s_ins = s_2;
        }
      } 
      else if( mode == -1 ) {  // else take the intersection with largest s 
        if( s_1 > s_2 ){
          //                                    streamlog_out(DEBUG0) << "ILDParallelPlanarMeasLayer::CalcXingPointWith: take 1st solution " << std::endl;
          x_ins = x_ins1;
          y_ins = y_ins1;
          s_ins = s_1;
        }
        else{
          //                                    streamlog_out(DEBUG0) << "ILDParallelPlanarMeasLayer::CalcXingPointWith: take 2nd solution " << std::endl;
          x_ins = x_ins2;
          y_ins = y_ins2;
          s_ins = s_2;
        }
      }
    }
  }
  


//  streamlog_out(DEBUG0) << "ILDParallelPlanarMeasLayer::CalcXingPointWith mode = " << mode << " dr = " << dr << " dz = " << z0 << " x = " << ref_point.x() << " y = "<< ref_point.y() << " z = " << ref_point.z() << " r = " << ref_point.Perp() << std::endl;
  
//  TVector3 xx_n;
//  int cuts = TVSurface::CalcXingPointWith(hel, xx_n, phi, mode, eps);
//  streamlog_out(DEBUG0) << "ILDParallelPlanarMeasLayer::CalcXingPointWith from Newton: cuts = " << cuts << " x = " << xx_n.x() << " y = "<< xx_n.y() << " z = " << xx_n.z() << " r = " << xx_n.Perp() << " phi = " << xx_n.Phi() << " dphi = " <<  phi << " " << this->TVMeasLayer::GetName() << std::endl;

  xx.SetXYZ(x_ins, y_ins, z_pca + s_ins * tanl);
  
  phi = -s_ins * omega ;
  
  // streamlog_out(DEBUG0) << "ILDParallelPlanarMeasLayer::CalcXingPointWith:             cuts = " << (IsOnSurface(xx) && (chg*phi*mode)<0)
  // << " x = " << xx.X()
  // << " y = " << xx.Y()
  // << " z = " << xx.Z()
  // << " r = " << xx.Perp()
  // << " phi = " << xx.Phi()
  // << " dphi = " <<  phi
  // << " " << this->TVMeasLayer::GetName() 
  // << std::endl;

  if( mode!=0 && fabs(phi)>1.e-10){ // (+1,-1) = (fwd,bwd)
    if( chg*phi*mode > 0){
      return 0;
    }
  }
    
  return (IsOnSurface(xx) ? 1 : 0);
  
}




