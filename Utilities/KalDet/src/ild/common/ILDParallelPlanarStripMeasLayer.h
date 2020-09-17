#ifndef __ILDParallelPlanarStripMeasLayer__
#define __ILDParallelPlanarStripMeasLayer__

/** ILDParallelPlanarStripMeasLayer: User defined KalTest measurement layer class 
 *
 * @author S.Aplin DESY
 */

//#include "TKalMatrix.h"
//#include "TVector3.h"
//#include "TVTrackHit.h"


#include "ILDParallelPlanarMeasLayer.h"

class ILDParallelPlanarStripMeasLayer : public ILDParallelPlanarMeasLayer {
  
public:
  
  /** Constructor Taking inner and outer materials, distance and phi of plane pca to origin, B-Field, Sorting policy, plane transverse witdth and offset of centre, longitudinal width, whether the layer is sensitive, Cell ID, and an optional name */
  ILDParallelPlanarStripMeasLayer(TMaterial &min,
                             TMaterial &mout,
                             Double_t   r,
                             Double_t   phi,
                             Double_t   Bz,
                             Double_t   SortingPolicy,
                             Double_t   xiwidth,
                             Double_t   zetawidth,
                             Double_t   xioffset,
                             Double_t   zoffset,
                             Double_t   UOrigin,
                             Double_t   stripAngle,
                             Int_t      CellID = -1,
                             const Char_t    *name = "ILDParallelPlanarStripMeasLayer")
  :
  ILDParallelPlanarMeasLayer(min,mout,r,phi,Bz,SortingPolicy,xiwidth,zetawidth,xioffset,zoffset,UOrigin,true,CellID,name), _stripAngle(stripAngle)
  
  { /* no op */ }
  
  
  // Parent's pure virtuals that must be implemented

  TKalMatrix XvToMv(const TVector3 &xv) const;

  TKalMatrix XvToMv(const TVTrackHit &, const TVector3   &xv) const {
    return XvToMv(xv);
  }

  TVector3 HitToXv(const TVTrackHit &vht) const ;
  
  void CalcDhDa(const TVTrackHit &vht, const TVector3   &xxv, const TKalMatrix &dxphiada, TKalMatrix &H)  const;
    
  ILDVTrackHit* ConvertLCIOTrkHit(edm4hep::ConstTrackerHit trkhit) const;
  
private:
  
  double _stripAngle;
  
};

#endif

