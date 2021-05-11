#ifndef ILDCYLINDERMEASLAYER_H
#define ILDCYLINDERMEASLAYER_H

/** ILDCylinderMeasLayer: User defined KalTest measurement layer class 
 *
 * @author S.Aplin DESY
 */


#include "ILDVMeasLayer.h"
#include <iostream>
#include <cmath>
/* #include "streamlog/streamlog.h" */


class ILDCylinderMeasLayer : public ILDVMeasLayer, public TCylinder {
  
public:
  
  /** Constructor Taking inner and outer materials, radius and half length, B-Field, whether the layer is sensitive, Cell ID, and an optional name */
  ILDCylinderMeasLayer(TMaterial &min,
                       TMaterial &mout,
                       Double_t   r0,
                       Double_t   lhalf,
                       Double_t   x0,
                       Double_t   y0,
                       Double_t   z0,
                       Double_t   Bz,
                       Bool_t     is_active,
                       Int_t      CellID = -1,
                       const Char_t    *name = "ILDCylinderMeasL") 
  : ILDVMeasLayer(min, mout, Bz, is_active, CellID, name),
  TCylinder(r0, lhalf,x0,y0,z0)
  { /* no op */ }
  

  Bool_t IsOnSurface(const TVector3 &xx) const {

    bool z = (xx.Z() >= GetZmin() && xx.Z() <= GetZmax());
    bool r = std::fabs( (xx-this->GetXc()).Perp() - this->GetR() ) < 1.e-3; // for very short, very stiff tracks this can be poorly defined, so we relax this here a bit to 1 micron

//    streamlog_out(DEBUG0) << "ILDCylinderMeasLayer IsOnSurface for " << this->TVMeasLayer::GetName() << " R =  " << this->GetR() << "  GetZmin() = " << GetZmin() << " GetZmax() = " << GetZmax()
//    << " dr = " << std::fabs( (xx-this->GetXc()).Perp() - this->GetR() ) << " r = " << r << " z = " << z 
//    << std::endl;
    
    return r && z;
  }

  

  // Parent's pure virtuals that must be implemented
  
  /** Global to Local coordinates */
  virtual TKalMatrix XvToMv    (const TVector3   &xv)   const;
  
  /** Global to Local coordinates */
  virtual TKalMatrix XvToMv    (const TVTrackHit &ht,
                                const TVector3   &xv)   const 
  
  { return this->XvToMv(xv); }  


  /** Local to Global coordinates */
  virtual TVector3   HitToXv   (const TVTrackHit &ht)   const;
  
  
  /** Calculate Projector Matrix */
  virtual void       CalcDhDa  (const TVTrackHit &ht,
                                const TVector3   &xv,
                                const TKalMatrix &dxphiada,
                                TKalMatrix &H)    const;
  
  /** Convert LCIO Tracker Hit to an ILDCylinderHit  */
  virtual ILDVTrackHit* ConvertLCIOTrkHit(edm4hep::ConstTrackerHit trkhit) const ;
  
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
  
};
#endif
