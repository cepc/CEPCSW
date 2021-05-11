#ifndef __ILDPLANARMEASLAYER__
#define __ILDPLANARMEASLAYER__
/** ILDRotatedTrapMeaslayer: User defined Rotated Trapezoid Planar KalTest measurement layer class used with ILDPLanarTrackHit.
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

class ILDRotatedTrapMeaslayer : public ILDVMeasLayer, public TPlane {
  
public:
  
  /** Constructor Taking inner and outer materials, centre and normal of the plane, B-Field, sorting policy, height, inner base length, outer base length, tilt angle around axis of symmetry, represents full or half petal, whether the layer is sensitive, Cell ID, and an optional name */
  ILDRotatedTrapMeaslayer(TMaterial &min,
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
                          Int_t      CellID = -1,
                          const Char_t    *name = "ILDRotatedTrapMeasL");
  
  
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
  
  /** Convert LCIO Tracker Hit to an ILDPLanarTrackHit  */
  virtual ILDVTrackHit* ConvertLCIOTrkHit(edm4hep::ConstTrackerHit trkhit) const ;
  
  /** Check if global point is on surface  */
  inline virtual Bool_t   IsOnSurface (const TVector3 &xx) const;
  
  /** Get sorting policy for this plane  */
  Double_t GetSortingPolicy() const { return _sortingPolicy; }
  
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
  Double_t _sortingPolicy ;
  Double_t _signZ ;
  Double_t _innerR ;
  Double_t _outerR ;
  Double_t _innerBaseLength ;
  Double_t _outerBaseLength ;
  Double_t _cosPhi ;  //** cos of the azimuthal angle of the petal 
  Double_t _sinPhi ;  //** sin of the azimuthal angle of the petal 
  Double_t _cosAlpha ; //** cos of the tilt angle of the petal 
  Double_t _sinAlpha ; //** sin of the tilt angle of the petal 
  Double_t _tanBeta ; //** tan of the openning angle of the petal
  
  // meaning of _halfPetal:
  //                  0 complete trapezoid
  //                 +1 positive half only, i.e. the side of the petal in which the transverse coordinate, mv(0,0), is positive 
  //                 -1 negative half only, i.e. the side of the petal in which the transverse coordinate, mv(0,0), is negative
  Int_t _halfPetal ;
  
  
  
};

#endif
