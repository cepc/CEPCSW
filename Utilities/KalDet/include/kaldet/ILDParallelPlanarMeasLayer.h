#ifndef __ILDParallelPlanarMeasLayer__
#define __ILDParallelPlanarMeasLayer__

/** ILDParallelPlanarMeasLayer: User defined KalTest measurement layer class 
 *
 * @author S.Aplin DESY
 */


#include "ILDPlanarMeasLayer.h"

class ILDParallelPlanarMeasLayer : public ILDPlanarMeasLayer {
  
public:
  
  /** Constructor Taking inner and outer materials, distance and phi of plane pca to origin, B-Field, Sorting policy, plane transverse witdth and offset of centre, longitudinal width, whether the layer is sensitive, Cell ID, and an optional name */
  ILDParallelPlanarMeasLayer(TMaterial &min,
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
                             Bool_t     is_active,
                             Int_t      CellID = -1,
                             const Char_t    *name = "ILDParallelPlanarMeasLayer")
  :
  ILDPlanarMeasLayer(min,mout,TVector3(r*cos(phi),r*sin(phi),zoffset),TVector3(cos(phi),sin(phi),0.0),Bz,SortingPolicy,xiwidth,zetawidth,xioffset,UOrigin,is_active,CellID,name), _r(r),_phi(phi),_cos_phi(cos(_phi)),_sin_phi(sin(_phi))
  { /* no op */ }
  
  
  // Parent's pure virtuals that must be implemented
  
  /** overloaded version of CalcXingPointWith using closed solution */
  virtual Int_t    CalcXingPointWith(const TVTrack  &hel,
                                     TVector3 &xx,
                                     Double_t &phi,
                                     Int_t     mode,
                                     Double_t  eps = 1.e-8) const;
  
  /** overloaded version of CalcXingPointWith using closed solution */
  virtual Int_t    CalcXingPointWith(const TVTrack  &hel,
                                     TVector3 &xx,
                                     Double_t &phi,
                                     Double_t  eps = 1.e-8) const{
    
    return CalcXingPointWith(hel,xx,phi,0,eps);
    
  }
  
  /** Get the intersection and the CellID, needed for multilayers */
  virtual int getIntersectionAndCellID(const TVTrack  &hel,
                                       TVector3 &xx,
                                       Double_t &phi,
                                       Int_t    &CellID,
                                       Int_t     mode,
                                       Double_t  eps = 1.e-8) const {
  
    CellID = this->getCellIDs()[0]; // not multilayer
    return CalcXingPointWith(hel,xx,phi,0,eps);
  
  }
  
  
protected:
  
  Double_t _r;
  Double_t _phi;
  Double_t _cos_phi;
  Double_t _sin_phi;
  
  
};

#endif

