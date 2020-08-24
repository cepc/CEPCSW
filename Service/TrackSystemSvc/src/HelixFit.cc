//
//  HelixFit.cc
//  MarlinTrk
//
//  Created by Steve Aplin on 9/16/11.
//

#include "TrackSystemSvc/IMarlinTrack.h"

#include "TrackSystemSvc/HelixFit.h"
#include "TrackSystemSvc/HelixTrack.h"
//#include "streamlog/streamlog.h"

#include "CLHEP/Matrix/SymMatrix.h"

namespace MarlinTrk{
  
  
  int HelixFit::fastHelixFit(int npt, double* xf, double* yf, float* rf, float* pf, double* wf, float* zf , float* wzf,int iopt, 
                             float* vv0, float* ee0, float& ch2ph, float& ch2z){
    
    
    
    if (npt < 3) {
      //streamlog_out(ERROR) << "Cannot fit less than 3 points return 1" << std::endl;
      ch2ph = 1.0e30;
      ch2z  = 1.0e30;
      return 1;
    }
    
    double eps = 1.0e-16;
#define ITMAX 15
#define MPT   600
#define MAX_CHI2 5000.0
    
    
    float sp2[MPT],
    del[MPT],deln[MPT],delzn[MPT],sxy[MPT],ss0[MPT],eee[MPT],
    delz[MPT],grad[5],cov[15],vv1[5],dv[5];
    
    //  double xf[MPT],yf[MPT],wf[MPT],zf[MPT],wzf[MPT];
    
    double alf,a0,a1,a2,a22,bet,cur,
    dd,den,det,dy,d2,f,fact,fg,f1,g,gam,gam0,g1,
    h,h2,p2,q2,rm,rn,
    xa,xb,xd,xi,xm,xx,xy,x1,x2,den2,
    ya,yb,yd,yi,ym,yy,y1,y2,wn,sa2b2,dd0,phic,aaa;
    
    
    xm = 0.;
    ym = 0.;
    
    for(int i=1; i < 15; ++i){
      ee0[i]=0.0;
    }
    
    for( int i=1; i < 5; ++i){
      grad[i]= 0.0;
      vv0[i] = 0.0;
    }
    
    float chi2=0.0;
    ch2ph = 0.0;
    ch2z = 0.0;
    
    for (int i = 0; i<npt; ++i) {
      sp2[i] = wf[i]*(rf[i]*rf[i]);
    }
    
    
    //std::cout <<  "npt = " <<  npt << std::endl;
    
    //  for(int i = 0; i<n; ++i){
    //std::cout << "xf = " <<  xf[i] << " " <<  i << std::endl; 
    //std::cout << "yf = " <<  yf[i] << " " <<  i << std::endl; 
    //std::cout << "zf = " <<  zf[i] << " " <<  i << std::endl; 
    //std::cout << "rf = " <<  rf[i] << " " <<  i << std::endl; 
    //std::cout << "pf = " <<  pf[i] << " " <<  i << std::endl; 
    //std::cout << "wf = " <<  wf[i] << " " <<  i << std::endl; 
    //std::cout << "wzf = " <<  wzf[i] << " " <<  i << std::endl; 
    //}
    
    
    wn=0.0;
    
    for (int i =0; i<npt; ++i) {
      xm = xm + xf[i]*wf[i];
      ym = ym + yf[i]*wf[i];
      wn = wn + wf[i];
    }
    
    rn = 1.0/wn;
    
    xm = xm * rn;
    ym = ym * rn;
    x2 = 0.;
    y2 = 0.;
    xy = 0.;
    xd = 0.;
    yd = 0.;
    d2 = 0.;
    
    //std::cout << "xm = " <<  xm <<  std::endl; 
    //std::cout << "ym = " <<  ym <<  std::endl; 
    //std::cout << "rn = " <<  rn <<  std::endl; 
    
    for (int i =0; i<npt; ++i) {
      xi = xf[i] - xm;
      yi = yf[i] - ym;
      xx = xi*xi;
      yy = yi*yi;
      x2 = x2 + xx*wf[i];
      y2 = y2 + yy*wf[i];
      xy = xy + xi*yi*wf[i];
      dd = xx + yy;
      xd = xd + xi*dd*wf[i];
      yd = yd + yi*dd*wf[i];
      d2 = d2 + dd*dd*wf[i];
    }
    
    
    //std::cout << "xd = " <<  xd <<  std::endl; 
    //std::cout << "yd = " <<  yd <<  std::endl; 
    //std::cout << "d2 = " <<  d2 <<  std::endl; 
    
    
    x2 = x2*rn;
    y2 = y2*rn;
    xy = xy*rn;
    d2 = d2*rn;
    xd = xd*rn;
    yd = yd*rn;
    f = 3.0* x2 + y2;
    g = 3.0* y2 + x2;
    fg = f*g;
    h = xy + xy;
    h2 = h*h;
    p2 = xd*xd;
    q2 = yd*yd;
    gam0 = x2 + y2;
    fact = gam0*gam0;
    a2 = (fg-h2-d2)/fact;
    fact = fact*gam0;
    a1 = (d2*(f+g) - 2.0*(p2+q2))/fact;
    fact = fact*gam0;
    a0 = (d2*(h2-fg) + 2.0*(p2*g + q2*f) - 4.0*xd*yd*h)/fact;
    a22 = a2 + a2;
    yb = 1.0e30;
    xa = 1.0;
    
    //std::cout << "a0 = " <<  a0 <<  std::endl; 
    //std::cout << "a22 = " <<  a22 <<  std::endl; 
    
    
    //  **                main iteration
    
    for (int i = 0 ; i < ITMAX; ++i) {
      
      ya = a0 + xa*(a1 + xa*(a2 + xa*(xa-4.0)));
      dy = a1 + xa*(a22 + xa*(4.0*xa - 12.0));
      xb = xa - ya/dy;
      
      if (fabs(ya) > fabs(yb)) {
        xb = 0.5 * (xb+xa) ;
      }
      
      if (fabs(xa-xb) < eps) break ;
      xa = xb;
      yb = ya;
      
      
    }
    
    //  **      
    
    gam = gam0*xb;
    f1 = f - gam;
    g1 = g - gam;
    x1 = xd*g1 - yd*h;
    y1 = yd*f1 - xd*h;
    det = f1*g1 - h2;
    den2= 1.0/(x1*x1 + y1*y1 + gam*det*det);
    
    if(den2 <= 0.0) {
      //                        streamlog_out(ERROR) << "den2 less than or equal to zero" 
      //                        << " x1 = " << x1 
      //                        << " y1 = " << y1
      //                        << " gam = " << gam  
      //                        << " det = " << det  
      //                        << std::endl;
      ch2ph = 1.0e30;
      ch2z  = 1.0e30;
      return 1;
    }
    
    den = sqrt(den2);
    cur = det*den + 0.0000000001 ;
    alf = -(xm*det + x1)*den ;
    bet = -(ym*det + y1)*den ;
    rm = xm*xm + ym*ym ;
    gam = ((rm-gam)*det + 2.0*(xm*x1 + ym*y1))*den*0.5;
    
    
    //  
    //  --------> calculation of standard circle parameters
    //            nb: cur is always positive
    
    double rr0=cur;
    double asym = bet*xm-alf*ym;
    double sst = 1.0;
    
    if(asym < 0.0) { 
      sst = -1.0 ;
    }
    
    rr0 = sst*cur;
    
    if( (alf*alf+bet*bet) <= 0.0 ){
      //                        streamlog_out(ERROR) << "(alf*alf+bet*bet) less than or equal to zero" << std::endl;
      ch2ph = 1.0e30;
      ch2z  = 1.0e30;
      return 1;
    }
    
    sa2b2 = 1.0/sqrt(alf*alf+bet*bet);
    dd0 = (1.0-1.0/sa2b2)/cur;
    aaa = alf*sa2b2;
    
    if( aaa >  1.0 ) aaa = 1.0;
    if( aaa < -1.0 ) aaa =-1.0;
    
    //  std::cout << std::setprecision(10) << "aaa = " <<  aaa <<  std::endl; 
    
    phic = asin(aaa)+ M_PI_2;
    
    if( bet > 0 ) phic = 2*M_PI - phic;
    
    double ph0 = phic + M_PI_2;
    
    if(rr0 <= 0.0)   ph0=ph0-M_PI;
    if(ph0 > 2*M_PI) ph0=ph0-2*M_PI;
    if(ph0 < 0.0)    ph0=ph0+2*M_PI;
    
    vv0[0] = rr0;
    vv0[2] = ph0;
    vv0[3] = dd0;
    
    //  std::cout << std::setprecision(10) << "rr0 = " <<  rr0 <<  std::endl; 
    //  std::cout << "ph0 = " <<  ph0 <<  std::endl; 
    //  std::cout << "dd0 = " <<  dd0 <<  std::endl; 
    
    double check=sst*rr0*dd0;
    //  std::cout << "check = " <<  check <<  std::endl; 
    
    if(check > 1.0-eps && check < 1.0+eps) {
      dd0 = dd0 - 0.007;
      vv0[3]=dd0;
    }
    
    //  
    //  -----> calculate phi distances to measured points
    //  
    
    double aa0 = sst;
    double ome = rr0;
    double gg0 = ome*dd0-aa0;
    
    double hh0 = 1.0/gg0;
    
    //  std::cout << "dd0 = " <<  dd0 <<  std::endl; 
    //  std::cout << "ome = " <<  ome <<  std::endl; 
    //  std::cout << "aa0 = " <<  aa0 <<  std::endl; 
    //  std::cout << "hh0 = " <<  hh0 <<  std::endl; 
    //  std::cout << "gg0 = " <<  gg0 <<  std::endl; 
    
    for (int i = 0 ; i < npt ; ++i) {
      
      asym   = bet*xf[i]-alf*yf[i] ;
      ss0[i] = 1.0 ;
      
      if(asym < 0.0) ss0[i] = -1.0 ;
      
      double ff0 = ome*(rf[i]*rf[i]-dd0*dd0)/(2.0*rf[i]*gg0) + dd0/rf[i];
      
      if(ff0 < -1.0) ff0 = -1.0;
      if(ff0 >  1.0) ff0 =  1.0;
      
      del[i]= ph0 + (ss0[i]-aa0)* M_PI_2 + ss0[i]*asin(ff0) - pf[i];
      
      //                std::cout << "asin(ff0) = " << asin(ff0)  << " i = " << i <<  std::endl; 
      //                std::cout << "aa0 = " <<  aa0 << " i = " << i <<  std::endl; 
      //                std::cout << "M_PI_2 = " <<  M_PI_2 << " i = " << i <<  std::endl; 
      //                std::cout << "pf[i] = " <<  pf[i] << " i = " << i <<  std::endl; 
      //                std::cout << "ff0 = " <<  ff0 << " i = " << i <<  std::endl; 
      //                std::cout << "ss0[i] = " <<  ss0[i] << " i = " << i <<  std::endl; 
      //                std::cout << "ph0 + (ss0[i]-aa0)* M_PI_2 = " <<  ph0 + (ss0[i]-aa0)* M_PI_2  << " i = " << i <<  std::endl; 
      //                std::cout << "ss0[i]*asin(ff0) = " <<  ss0[i]*asin(ff0)  << " i = " << i <<  std::endl; 
      //                std::cout << "ph0 + (ss0[i]-aa0)* M_PI_2 + ss0[i]*asin(ff0) = " << ph0 + (ss0[i]-aa0)* M_PI_2 + ss0[i]*asin(ff0) << std::endl;
      //                std::cout << "del[i] = " <<  del[i] << " i = " << i <<  std::endl; 
      
      if(del[i] >  M_PI) del[i] = del[i] - 2*M_PI;
      if(del[i] < -M_PI) del[i] = del[i] + 2*M_PI;
      
      
      
    }
    
    
    
    //  
    //  -----> fit straight line in s-z
    //  
    
    
    for (int i = 0 ; i < npt ; ++i) {
      
      eee[i] = 0.5*vv0[0] * sqrt( fabs( (rf[i] * rf[i]-vv0[3]*vv0[3]) / (1.0-aa0*vv0[0]*vv0[3]) ) );
      
      if(eee[i] >  0.99990)  eee[i]=  0.99990;
      if(eee[i] < -0.99990)  eee[i]= -0.99990;
      
      sxy[i]=2.0*asin(eee[i])/ome;
      
    }
    
    
    double sums  = 0.0;
    double sumss = 0.0;
    double sumz  = 0.0;
    double sumzz = 0.0;
    double sumsz = 0.0;
    double sumw  = 0.0;
    
    for (int i = 0; i<npt; ++i) {
      sumw  = sumw  +                 wzf[i];
      sums  = sums  + sxy[i]        * wzf[i];
      sumss = sumss + sxy[i]*sxy[i] * wzf[i];
      sumz  = sumz  + zf[i]         * wzf[i];
      sumzz = sumzz + zf[i]*zf[i]   * wzf[i];
      sumsz = sumsz + zf[i]*sxy[i]  * wzf[i];
    }
    
    double denom = sumw*sumss - sums*sums;
    
    if (fabs(denom) < eps){
      //                        streamlog_out(ERROR) << "fabs(denom) less than or equal to zero" << std::endl;
      ch2ph = 1.0e30;
      ch2z  = 1.0e30;
      return 1;
    }
    
    double dzds  = (sumw*sumsz-sums*sumz) /denom;
    double zz0   = (sumss*sumz-sums*sumsz)/denom;
    
    vv0[1]= dzds;
    vv0[4]= zz0;
    
    //  
    //  -----> calculation chi**2
    //  
    for (int i = 0 ; i<npt; ++i) {
      
      delz[i]= zz0+dzds*sxy[i]-zf[i];
      ch2ph = ch2ph + sp2[i]*del[i]*del[i];
      ch2z = ch2z + wzf[i]*delz[i]*delz[i];
      chi2 = ch2ph + ch2z;
      
    }
    
    if(chi2 > MAX_CHI2) {
      //                        streamlog_out(ERROR) << "Chi2 greater than " <<  MAX_CHI2 << "return 1 " << std::endl;
      ch2ph = 1.0e30;
      ch2z  = 1.0e30;
      return 1;
    }
    
    for (int i = 0 ; i<npt; ++i) {
      
      double ff0 = ome*(rf[i]*rf[i]-dd0*dd0)/(2.0*rf[i]*gg0) + dd0/rf[i];
      
      if (ff0 >  0.99990) ff0 =  0.99990;
      
      if (ff0 < -0.99990) ff0 = -0.99990;
      
      double eta = ss0[i]/sqrt(fabs((1.0+ff0)*(1.0-ff0)));
      double dfd = (1.0+hh0*hh0*(1.0-ome*ome*rf[i]*rf[i]))/(2.0*rf[i]);
      double dfo = -aa0*(rf[i]*rf[i]-dd0*dd0)*hh0*hh0/(2.0*rf[i]);
      double dpd = eta*dfd;
      double dpo = eta*dfo;
      
      //        -----> derivatives of z component
      double ggg = eee[i]/sqrt(fabs( (1.0+eee[i])*(1.0-eee[i])));
      double dza = sxy[i];
      check = rf[i]*rf[i]-vv0[3]*vv0[3];
      
      if(fabs(check) > 1.0-eps) check=2.*0.007;
      
      double dzd = 2.0*( vv0[1]/vv0[0] ) * fabs( ggg ) * ( 0.5*aa0*vv0[0]/( 1.0-aa0*vv0[3]*vv0[0] )-vv0[3]/check );
      
      double dzo = -vv0[1]*sxy[i]/vv0[0] + vv0[1] * ggg/( vv0[0]*vv0[0]) * ( 2.0+ aa0*vv0[0]*vv0[3]/(1.0-aa0*vv0[0]*vv0[3]) );
      
      //  -----> error martix
      
      ee0[0] = ee0[0] + sp2[i]*  dpo*dpo  + wzf[i] * dzo*dzo;
      ee0[1] = ee0[1]                     + wzf[i] * dza*dzo;
      ee0[2] = ee0[2]                     + wzf[i] * dza*dza;
      ee0[3] = ee0[3] + sp2[i]*  dpo;
      ee0[4] = ee0[4];
      ee0[5] = ee0[5] + sp2[i];
      ee0[6] = ee0[6] + sp2[i]*  dpo*dpd  + wzf[i] * dzo*dzd;
      ee0[7] = ee0[7]                     + wzf[i] * dza*dzd;
      ee0[8] = ee0[8] + sp2[i]*      dpd;
      ee0[9] = ee0[9] + sp2[i]*  dpd*dpd  + wzf[i] * dzd*dzd;
      ee0[10]= ee0[10]                    + wzf[i] * dzo;
      ee0[11]= ee0[11]                    + wzf[i] * dza;
      ee0[12]= ee0[12];
      ee0[13]= ee0[13]                    + wzf[i] * dzd;
      ee0[14]= ee0[14]                    + wzf[i];
      
      //        -----> gradient vector
      grad[0]=grad[0] - del[i] *sp2[i]*dpo - delz[i]*wzf[i]*dzo;
      grad[1]=grad[1] -                      delz[i]*wzf[i]*dza;
      grad[2]=grad[2] - del[i] *sp2[i];
      grad[3]=grad[3] - del[i] *sp2[i]*dpd - delz[i]*wzf[i]*dzd;
      grad[4]=grad[4] -                      delz[i]*wzf[i];
      
      
    }
    
    if (iopt < 3) {
      /*
      streamlog_out(DEBUG1) << "HelixFit: " <<
      " d0 = " << vv0[3] <<
      " phi0 = " << vv0[2] <<
      " omega = " << vv0[0] <<
      " z0 = " << vv0[4] <<
      " tanL = " << vv0[1] <<
      std::endl;
      */
      return 0 ;
    }
    
    // ------> NEWTONS NEXT GUESS       
    
    for (int i =0; i<15; ++i) {
      cov[i]=ee0[i];
    }
    
    // Convert Covariance Matrix
    CLHEP::HepSymMatrix cov0(5) ; 
    
    int icov = 0 ;
    
    for(int irow=0; irow<5; ++irow ){
      for(int jcol=0; jcol<irow+1; ++jcol){
        //      std::cout << "row = " << irow << " col = " << jcol << std::endl ;
        //      std::cout << "cov["<< icov << "] = " << _cov[icov] << std::endl ;
        cov0[irow][jcol] = ee0[icov] ;
        ++icov ;
      }
    }
    
    int error = 0 ;
    cov0.invert(error);
    
    if( error != 0 ) {
      //streamlog_out(ERROR) << "CLHEP Matrix inversion failed" << "return 1 " << std::endl;
      ch2ph = 1.0e30;
      ch2z  = 1.0e30;
      return 1;
    }
    
    icov = 0 ;
    
    for(int irow=0; irow<5; ++irow ){
      for(int jcol=0; jcol<irow+1; ++jcol){
        //      std::cout << "row = " << irow << " col = " << jcol << std::endl ;
        //      std::cout << "cov["<< icov << "] = " << _cov[icov] << std::endl ;
        cov[icov] = cov0[irow][jcol];
        ++icov ;
      }
    }
    
    // here we need to invert the matrix cov[15]
    
    
    
    for (int i = 0; i<5; ++i) {
      dv[i] = 0.0 ;
      for (int j = 0; j<5; ++j) {
        int index = 0;
        if ( i >= j ) {
          index = (i*i-i)/2+j ;
        }
        else{
          index = (j*j-j)/2+i;
        }
        dv[i] = dv[i] + cov[index] * grad[j]    ;
      }
    }
    
    
    for (int i = 0; i<5; ++i) {
      vv1[i] = vv0[i]+dv[i];
    }
    
    gg0 = vv1[0]*vv1[3]-aa0;
    
    for (int i=0; i<npt; ++i) {
      double ff0 = vv1[0]*(rf[i]*rf[i]-vv1[3]*vv1[3]) / (2.0*rf[i]*gg0) + vv1[3]/rf[i];
      
      if (ff0 >  1) ff0 =  1.0;
      if (ff0 < -1) ff0 = -1.0;
      
      deln[i] = vv1[2] + (ss0[i]-aa0)*M_PI_2+ss0[i]*asin(ff0)-pf[i];
      
      if(deln[i] >  M_PI) deln[i] = deln[i] - 2*M_PI;
      if(deln[i] < -M_PI) deln[i] = deln[i] + 2*M_PI;
      
      eee[i] = 0.5*vv1[0]*sqrt( fabs( (rf[i]*rf[i]-vv1[3]*vv1[3]) / (1.0-aa0*vv1[0]*vv1[3]) ) );
      
      if(eee[i] >  0.99990)  eee[i]=  0.99990;
      if(eee[i] < -0.99990)  eee[i]= -0.99990;
      
      sxy[i] = 2.0*asin(eee[i])/vv1[0];
      delzn[i]= vv1[4]+vv1[1]*sxy[i]-zf[i];
      
    }
    
    double chi1 = 0.0;
    ch2ph = 0.0;
    ch2z = 0.0;
    
    for (int i =0; i<5; ++i) {
      chi1   = chi1  + sp2[i]*deln[i]*deln[i] + wzf[i]*delzn[i]*delzn[i];
      ch2ph  = ch2ph + sp2[i]*deln[i]*deln[i];
      ch2z   = ch2z  + wzf[i]*delzn[i]*delzn[i];        
    }
    
    if (chi1<chi2) {
      for (int i =0; i<5; ++i) {
        vv0[i] = vv1[i];
      }
      chi2 = chi1;      
    }
    std::cout << "HelixFit: " <<    
      //streamlog_out(DEBUG1) << "HelixFit: " <<
      " d0 = " << vv0[3] <<
      " phi0 = " << vv0[2] <<
      " omega = " << vv0[0] <<
      " z0 = " << vv0[4] <<
      " tanL = " << vv0[1] <<
      std::endl;
                                    
    return 0;
    
  }
  
  
}

















