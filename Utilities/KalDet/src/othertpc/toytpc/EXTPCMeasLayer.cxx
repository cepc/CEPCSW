//*************************************************************************
//* ===================
//*  EXTPCMeasLayer Class
//* ===================
//*
//* (Description)
//*   Sample measurement layer class used by EXTPCHit.
//* (Requires)
//* (Provides)
//*     class EXTPCMeasLayer
//* (Update Recored)
//*   2003/09/30  Y.Nakashima       Original version.
//*
//*************************************************************************
//

#include "EXTPCMeasLayer.h"
#include "EXTPCHit.h"
#include "EXTPCKalDetector.h"
#include "TRandom.h"
#include "TMath.h"


ClassImp(EXTPCMeasLayer)
                                                                                
EXTPCMeasLayer::EXTPCMeasLayer(TMaterial &min,
                               TMaterial &mout,
                	       Double_t   y0,
                               Double_t   lhalf,
                               Double_t   sigmax0,
                               Double_t   sigmax1,
                               Double_t   sigmaz0,
                               Double_t   sigmaz1,
                               Bool_t     type)
              : EXVMeasLayer(min, mout, type),
		TPlane(TVector3(0.,y0,lhalf), TVector3(0.,1.,0.)),
                fSigmaX0(sigmax0),
                fSigmaX1(sigmax1),
                fSigmaZ0(sigmaz0),
                fSigmaZ1(sigmaz1),
                fLayer(-1)
{
}

EXTPCMeasLayer::EXTPCMeasLayer(TMaterial &min,
                               TMaterial &mout,
			       Double_t   y0,
			       Double_t   XMin,
			       Double_t   XMax,
                               Double_t   lhalf,
                               Double_t   sigmax0,
                               Double_t   sigmax1,
                               Double_t   sigmaz0,
                               Double_t   sigmaz1,
                               Bool_t     type,
                               Int_t      layer,
                               Int_t      module)
              : EXVMeasLayer(min, mout, type),
		TPlane(TVector3(0.,y0,lhalf), TVector3(0.,1.,0.)),
		fXMin(XMin),
		fXMax(XMax),
                fSigmaX0(sigmax0),
                fSigmaX1(sigmax1),
                fSigmaZ0(sigmaz0),
                fSigmaZ1(sigmaz1),
                fModule(module),
                fLayer(layer)
{
}

EXTPCMeasLayer::~EXTPCMeasLayer()
{
}

TKalMatrix EXTPCMeasLayer::XvToMv(const TVector3 &xv,
                                        Int_t     side) const
{
  TKalMatrix mv(kMdim,1);
  mv(0,0) =xv.X();
  mv(1,0) =xv.Z();
  return mv;
}

TKalMatrix EXTPCMeasLayer::XvToMv(const TVTrackHit &vht,
                                  const TVector3   &xv) const
{
   return XvToMv(xv, dynamic_cast<const EXTPCHit &>(vht).GetSide());
}

TVector3 EXTPCMeasLayer::HitToXv(const TVTrackHit &vht) const
{
   const EXTPCHit &ht = dynamic_cast<const EXTPCHit &>(vht);
   Double_t x = ht(0,0);
   Double_t y = GetXc().Y();
   Double_t z = ht(1,0);
   return TVector3(x,y,z);
}

void EXTPCMeasLayer::CalcDhDa(const TVTrackHit &vht,
			      const TVector3 &xxv,
			      const TKalMatrix &dxphiada,
			      TKalMatrix &H) const
{
  Int_t sdim = H.GetNcols();
  Int_t hdim = TMath::Max(5,sdim-1);

  for (Int_t i=0;i<hdim;i++){
    H(0,i) = dxphiada(0,i);
    H(1,i) = dxphiada(2,i);
  }
}

/*
 void EXTPCMeasLayer::CalcDhDa(const TVTrackHit &vht,
                              const TVector3   &xxv,
                              const TKalMatrix &dxphiada,
                                    TKalMatrix &H)  const
 {
  const EXTPCHit &ht = dynamic_cast<const EXTPCHit &>(vht);

  //  Calculate
  //     H = (@h/@a) = (@phi/@a, @z/@a)^t
  //  where
  //         h(a) = (phi, z)^t: expected meas vector
  //         a = (drho, phi0, kappa, dz, tanl, t0)
 

   Int_t sdim = H.GetNcols();
      Int_t hdim = TMath::Max(5,sdim-1);

   TVector3 xxvc = xxv - GetXc();
   Double_t xv = xxvc.X();
   Double_t yv = xxvc.Y();
   Double_t xxyy = xv * xv + yv * yv;
   
   // Set H = (@h/@a) = (@d/@a, @z/@a)^t
   
   for (Int_t i=0; i<hdim; i++) {
      H(0,i) = - (yv / xxyy) * dxphiada(0,i) 
               + (xv / xxyy) * dxphiada(1,i);
      H(0,i) *= 0;//GetR();
#if 0
      H(1,i) = - ht.GetSide() *  dxphiada(2,i);
#else
      H(1,i) = dxphiada(2,i);
#endif
   }
   if (sdim == 6) {
      H(0,sdim-1) = 0.;
#if 0
      H(1,sdim-1) = ht.GetVdrift();
#else
      H(1,sdim-1) = - ht.GetVdrift();
#endif
   }
}
*/

Double_t EXTPCMeasLayer::GetSigmaX(Double_t zdrift) const
{
   return TMath::Sqrt(fSigmaX0 * fSigmaX0 + fSigmaX1 * fSigmaX1 * zdrift);
}

Double_t EXTPCMeasLayer::GetSigmaZ(Double_t zdrift) const
{
   return TMath::Sqrt(fSigmaZ0 * fSigmaZ0 + fSigmaZ1 * fSigmaZ1 * zdrift);
}

Double_t EXTPCMeasLayer::GetSortingPolicy() const
{
  return GetXc().Y();
}


void EXTPCMeasLayer::ProcessHit(const TVector3  &xx,
                                      TObjArray &hits)
{
   Int_t      side = (xx.Z() < 0. ? -1 : 1);
   TKalMatrix h    = XvToMv(xx, side);
   Double_t   rphi = h(0, 0);
   Double_t   d    = h(1, 0);

   Double_t dx = GetSigmaX(d);
   Double_t dz = GetSigmaZ(d);
   rphi += gRandom->Gaus(0., dx);   // smearing rphi
   d    += gRandom->Gaus(0., dz);   // smearing drift distance

   Double_t v = EXTPCKalDetector::GetVdrift();
#if 0
   d         += v * EXEventGen::GetT0(); // T0 shift
#endif

   Double_t meas [2];
   Double_t dmeas[2];
   meas [0] = rphi;
   meas [1] = d;
   dmeas[0] = dx;
   dmeas[1] = dz;

   Double_t b = EXTPCKalDetector::GetBfield();
   hits.Add(new EXTPCHit(*this, meas, dmeas, side, v, xx, b));
}
