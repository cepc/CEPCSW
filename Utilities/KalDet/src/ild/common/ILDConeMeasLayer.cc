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

#include "ILDConeMeasLayer.h"
#include "TMath.h"

 
ILDConeMeasLayer::ILDConeMeasLayer(TMaterial &min,
                                     TMaterial &mout,
                                     Double_t   z1,
                                     Double_t   r1,
                                     Double_t   z2,
                                     Double_t   r2,
                                     Double_t   Bz,
                                     Double_t   SortingPolicy,
                                     Bool_t     is_active,
                                     Int_t      CellID,
                               const Char_t    *name)
                   : ILDVMeasLayer(min, mout, Bz, is_active, CellID, name),
                   TCutCone(r1*(z2-z1)/(r2-r1), 
                            r2*(z2-z1)/(r2-r1), 
                               (r2-r1)/(z2-z1),
                            0.,0.,(r2*z1-r1*z2)/(r2-r1)),
                   fZ1(z1),
                   fR1(r1),
                   fZ2(z2),
                   fR2(r2),
                   fsortingPolicy(SortingPolicy)
{
}

ILDConeMeasLayer::~ILDConeMeasLayer()
{
}

TKalMatrix ILDConeMeasLayer::XvToMv(const TVector3 &xxv) const
{
   // Calculate hit coordinate information:
   //	mv(0,0) = r * phi 
   //     (1,0) = z

   TKalMatrix mv(kMdim,1);
   TVector3 xv = xxv - GetXc();
   Double_t r  = xv.Z()*GetTanA();
    
   mv(0,0)  = r * TMath::ATan2(xv.Y(), xv.X());
   mv(1,0)  = xv.Z();
   return mv;
}

TKalMatrix ILDConeMeasLayer::XvToMv(const TVTrackHit &,
                                     const TVector3   &xv) const
{
   return XvToMv(xv);
}

TVector3 ILDConeMeasLayer::HitToXv(const TVTrackHit &vht) const
{
//    const EXBPConeHit &ht = dynamic_cast<const EXBPConeHit &>(vht);
// 
//    Double_t r   = ht(1,0) * GetTanA();
//    Double_t phi = ht(0,0) / r;
//    Double_t x   = GetXc().X() + r * TMath::Cos(phi);
//    Double_t y   = GetXc().Y() + r * TMath::Sin(phi);
//    Double_t z   = GetXc().Z() + ht(1,0);

   // streamlog_out( ERROR ) << "Don't use this, it's not implemented!";

   return TVector3(0.,0.,0.);
}

void ILDConeMeasLayer::CalcDhDa(const TVTrackHit &vht,
                                 const TVector3   &xxv,
                                 const TKalMatrix &dxphiada,
                                       TKalMatrix &H)  const
{
   // Calculate
   //    H = (@h/@a) = (@phi/@a, @z/@a)^t
   // where
   //        h(a) = (phi, z)^t: expected meas vector
   //        a = (drho, phi0, kappa, dz, tanl, t0)
   //

   Int_t sdim = H.GetNcols();
   Int_t hdim = TMath::Max(5,sdim-1);

   TVector3 xv = xxv - GetXc();
   Double_t x  = xv.X();
   Double_t y  = xv.Y();
   Double_t z  = xv.Z();
   Double_t xxyy = x * x + y * y;
   Double_t phi  = TMath::ATan2(y, x);
   Double_t tana = GetTanA();
   Double_t r    = z * GetTanA();

   // Set H = (@h/@a) = (@d/@a, @z/@a)^t
   
   for (Int_t i=0; i<hdim; i++) {
      H(0,i) = - r * (y / xxyy) * dxphiada(0,i) 
               + r * (x / xxyy) * dxphiada(1,i)
               +     tana * phi * dxphiada(2,i);
      H(1,i) =  dxphiada(2,i);
   }
   if (sdim == 6) {
      H(0,sdim-1) = 0.;
      H(1,sdim-1) = 0.;
   }
}

Bool_t ILDConeMeasLayer::IsOnSurface(const TVector3 &xx) const
{
    TVector3 xxc = xx - GetXc();
    Double_t r   = xxc.Perp();
    Double_t z   = xxc.Z();
    Double_t s   = (r - GetTanA()*z) * (r + GetTanA()*z);
    const Double_t kTol = 1.e-8;

#if 0
    std::cout << this->TVMeasLayer::GetName() << ":" << this->GetIndex() << ":" << std::endl;
    std::cout << "s=" << s << " xx=(" << xx.X() << "," << xx.Y() << "," << xx.Z() << ")" << std::endl;
    std::cout << "bool=" << (TMath::Abs(s) < kTol && ((xx.Z()-fZ1)*(xx.Z()-fZ2) <= 0.)) << std::endl;
    std::cout << "fZ1=" << fZ1 << " fZ2=" << fZ2 << std::endl;
#endif

    return (TMath::Abs(s) < kTol && ((xx.Z()-fZ1)*(xx.Z()-fZ2) <= 0.));
} 

