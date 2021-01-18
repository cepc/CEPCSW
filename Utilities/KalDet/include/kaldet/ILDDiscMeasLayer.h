#ifndef __ILDDISCMEASLAYER__
#define __ILDDISCMEASLAYER__

/** ILDDiscMeasLayer: User defined KalTest Disc measurement layer class used with ILDPLanarTrackHit. 
 *
 * @author S.Aplin DESY
 */

#include "TVector3.h"

#include "kaltest/TKalMatrix.h"
#include "kaltest/TPlane.h"
#include "kaltest/KalTrackDim.h"
#include "ILDVMeasLayer.h"

#include "TMath.h"
#include <sstream>

class TVTrackHit;


class ILDDiscMeasLayer : public ILDVMeasLayer, public TPlane {
  
public:
  
  /** Constructor Taking inner and outer materials, center and normal to the plane, B-Field, sorting policy, min and max r, whether the layer is sensitive, Cell ID, and an optional name */
  
  ILDDiscMeasLayer(TMaterial &min,
                   TMaterial &mout,
                   const TVector3  &center,
                   const TVector3  &normal,
                   double   Bz,
                   double   SortingPolicy,
                   double   rMin,
                   double   rMax,
                   Bool_t     is_active,
                   Int_t      CellID = -1,
                   const Char_t    *name = "ILDDiscMeasL")
  : ILDVMeasLayer(min, mout, Bz, is_active, CellID, name),
  TPlane(center, normal),
  _sortingPolicy(SortingPolicy), _rMin(rMin), _rMax(rMax)
  { /* no op */ }
  
  
  
  // Parrent's pure virtuals that must be implemented
  
  /** Global to Local coordinates */
  virtual TKalMatrix XvToMv    (const TVTrackHit &ht,
                                const TVector3   &xv) const 
  { return this->XvToMv(xv); }
  
  /** Global to Local coordinates */
  virtual TKalMatrix XvToMv    (const TVector3   &xv) const;
  
  /** Local to Global coordinates */  
  virtual TVector3   HitToXv   (const TVTrackHit &ht) const;
  
  /** Calculate Projector Matrix */
  virtual void       CalcDhDa  (const TVTrackHit &ht,
                                const TVector3   &xv,
                                const TKalMatrix &dxphiada,
                                TKalMatrix &H)  const;
  
  virtual Int_t CalcXingPointWith(const TVTrack  &hel,
                                  TVector3 &xx,
                                  Double_t &phi,
                                  Int_t     mode,
                                  Double_t  eps) const;
    
  
  /** Convert LCIO Tracker Hit to an ILDPLanarTrackHit  */
  virtual ILDVTrackHit* ConvertLCIOTrkHit(edm4hep::ConstTrackerHit trkhit) const ;
  
  /** Check if global point is on surface  */
  inline virtual Bool_t   IsOnSurface (const TVector3 &xx) const;
  
  /** Get sorting policy for this plane  */
  double GetSortingPolicy() const { return _sortingPolicy; }
  
  /** Get the intersection and the CellID, needed for multilayers */
  virtual int getIntersectionAndCellID(const TVTrack  &hel,
                                       TVector3 &xx,
                                       Double_t &phi,
                                       Int_t    &CellID,
                                       Int_t     mode,
                                       Double_t  eps = 1.e-8) const {
    
    CellID = this->getCellIDs()[0]; // not multilayer
    return this->CalcXingPointWith(hel,xx,phi,0,eps);
  
  }
  
private:
  double _sortingPolicy;
  double _rMin;
  double _rMax;
  
};

#endif
