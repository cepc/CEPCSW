#include "Tools/FTDHelixFitter.h"

#include <sstream>
#include <algorithm>
#include <cmath>

#include "edm4hep/TrackerHit.h"
#include "UTIL/BitSet32.h"
#include "UTIL/ILDConf.h"
#include "UTIL/LCTrackerConf.h"
//#include "marlin/VerbosityLevels.h"
#include "TrackSystemSvc/HelixFit.h"

#include "Tools/KiTrackMarlinTools.h"


FTDHelixFitter::FTDHelixFitter( std::vector<edm4hep::ConstTrackerHit> trackerHits ){
  _trackerHits = trackerHits;
  fit();
}

FTDHelixFitter::FTDHelixFitter( edm4hep::Track* track ){
  _trackerHits.clear();
  //int nHits = track->trackerHits_size();
  std::copy(track->trackerHits_begin(), track->trackerHits_end(), std::back_inserter(_trackerHits));
  //for(int i=0;i<nHits;i++){
  //  edm4hep::ConstTrackerHit hit = &track->getTrackerHits(i);
  //  _trackerHits.push_back(hit);
  //}
  fit();
}

void FTDHelixFitter::fit(){
  std::sort( _trackerHits.begin(), _trackerHits.end(), KiTrackMarlin::compare_TrackerHit_z );
   
  int nHits = _trackerHits.size();
  int iopt = 2;
  float chi2RPhi;
  float chi2Z;
  
  if( nHits < 3 ){
    std::stringstream s;
    s << "FTDHelixFitter::fit(): Cannot fit less with less than 3 hits. Number of hits =  " << nHits << "\n";
      
    throw FTDHelixFitterException( s.str() );
  }
   
  double* xh  = new double[nHits];
  double* yh  = new double[nHits];
  float*  zh  = new float[nHits];
  double* wrh = new double[nHits];
  float*  wzh = new float[nHits];
  float*  rh  = new float[nHits];
  float*  ph  = new float[nHits];
  
  float par[5];
  float epar[15];
  
  for( int i=0; i<nHits; i++ ){
    edm4hep::ConstTrackerHit hit = _trackerHits[i];
      
    xh[i] = hit.getPosition()[0];
    yh[i] = hit.getPosition()[1];
    zh[i] = hit.getPosition()[2];
      
    float resZ = 0.1;
    wzh[i] = 1.0/( resZ * resZ );
      
    rh[i] = float(sqrt(xh[i]*xh[i]+yh[i]*yh[i]));
    ph[i] = atan2(yh[i],xh[i]);
    if (ph[i] < 0.)  ph[i] = 2.*M_PI + ph[i]; 
            
    if( UTIL::BitSet32( hit.getType() )[ UTIL::ILDTrkHitTypeBit::COMPOSITE_SPACEPOINT ] ){
      float sigX = hit.getCovMatrix()[0];
      float sigY = hit.getCovMatrix()[2];
      wrh[i] = 1/sqrt( sigX*sigX + sigY*sigY );
    }
    else {
      //TrackerHitPlane* hitPlane = dynamic_cast<TrackerHitPlane*>( hit );
      //wrh[i] = double(1.0/( hitPlane->getdU()*hitPlane->getdU() + hitPlane->getdV()*hitPlane->getdV() ) );
      wrh[i] = double(1.0/( hit.getCovMatrix(2)*hit.getCovMatrix(2) + hit.getCovMatrix(5)*hit.getCovMatrix(5) ));
    }
  }
     
  MarlinTrk::HelixFit helixFitter;
   
  helixFitter.fastHelixFit(nHits, xh, yh, rh, ph, wrh, zh, wzh,iopt, par, epar, chi2RPhi, chi2Z);
  par[3] = par[3]*par[0]/fabs(par[0]);
   
  _omega = par[0];
  _tanLambda = par[1];
  _phi0 = par[2];
  _d0 = par[3];
  _z0 = par[4];
  
  float chi2 = chi2RPhi+chi2Z;
  int Ndf = 2*nHits-5;
     
  delete[] xh;
  delete[] yh;
  delete[] zh;
  delete[] wrh;
  delete[] wzh;
  delete[] rh;
  delete[] ph;
  xh  = NULL; 
  yh  = NULL; 
  zh  = NULL; 
  wrh = NULL; 
  wzh = NULL; 
  rh = NULL; 
  ph = NULL; 
     
  //streamlog_out(DEBUG1) << "chi2 = " << chi2 << ", Ndf = " << Ndf << "\n";
  
  _chi2 = chi2;
  _Ndf = Ndf;
  
  return;
}




