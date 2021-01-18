#ifndef ILDCONEMEASLAYER_H
#define ILDCONEMEASLAYER_H
//*************************************************************************
//* ===================
//*  ILDConeMeasLayer Class
//* ===================
//*
//* (Update Recored)
//*   2012/01/19  K.Fujii       Original version. (EXBPConeMeasLayer)
//*   2012/01/24 R.Glattauer    Adapted to ILD common in KalDet
//*
//*************************************************************************
//
#include "TVector3.h"

#include "kaltest/TKalMatrix.h"
#include "kaltest/TCutCone.h"
#include "kaltest/KalTrackDim.h"

#include "ILDVMeasLayer.h"
#include "iostream"
/* #include "streamlog/streamlog.h" */
#include "UTIL/ILDConf.h"
#include "edm4hep/TrackerHit.h"

class ILDConeMeasLayer : public ILDVMeasLayer, public TCutCone {
public:
   // Ctors and Dtor
   /** Constructor Taking inner and outer materials, z and radius at start and end, B-Field, whether the layer is sensitive, Cell ID, and an optional name */
   ILDConeMeasLayer(TMaterial &min,
                     TMaterial &mout,
                     Double_t   z1,
                     Double_t   r1,
                     Double_t   z2,
                     Double_t   r2,
                     Double_t   Bz,
                     Double_t   SortingPolicy,
                     Bool_t     is_active,
                     Int_t      CellID = -1,
               const Char_t    *name = "BPCONEML");
   virtual ~ILDConeMeasLayer();

   // Parrent's pure virtuals that must be implemented
   /** Global to Local coordinates */
   virtual TKalMatrix XvToMv    (const TVTrackHit &ht,
                                 const TVector3   &xv) const;
   
   /** Global to Local coordinates */
   virtual TKalMatrix XvToMv    (const TVector3   &xv) const;
   
   /** Local to Global coordinates */
   virtual TVector3   HitToXv   (const TVTrackHit &ht) const;
   
   /** Calculate Projector Matrix */
   virtual void       CalcDhDa  (const TVTrackHit &ht,
                                 const TVector3   &xv,
                                 const TKalMatrix &dxphiada,
                                       TKalMatrix &H)  const;

   Bool_t IsOnSurface(const TVector3 &xx) const;

   /** Convert LCIO Tracker Hit to an ILDCylinderHit  */
   virtual ILDVTrackHit* ConvertLCIOTrkHit(edm4hep::ConstTrackerHit trkhit) const {
      
      /* streamlog_out( ERROR ) << "Don't use this, it's not implemented!"; */
      return NULL;
   }
   
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

  /** Get sorting policy for this plane  */
  virtual double GetSortingPolicy() const { return fsortingPolicy; }


  
   
private:
  Double_t fZ1;      // z of front face
  Double_t fR1;      // r of front face
  Double_t fZ2;      // z of back end
  Double_t fR2;      // r of back end
  Double_t fsortingPolicy; // used for sorting the layers in to out

  
};

#endif
