//*************************************************************************
//* ====================
//*  TVMeasLayer Class
//* ====================
//*
//* (Description)
//*   Measurement layer interface class.
//* (Requires)
//* (Provides)
//*     class TVMeasLayer
//* (Update Recored)
//*   2003/09/30  Y.Nakashima       Original version.
//*   2005/02/23  A.Yamaguchi       Added new data members, fFwdX0Inv,
//*                                 fBwdX0Inv and fIndex, and their
//*                                 corresponding getters and setters.
//*                                 Added a new method, GetX0Inv().
//*   2005/0X/XX  A.Yamaguchi       Replaced fFwdX0Inv, fBwdX0Inv, and
//*                                 their getters and setters by
//*                                 fMaterialOutPtr, fMaterialInPtr, and
//*                                 their getters and setters.
//*   2005/08/15  K.Fujii           Added fIsActive and IsActive().
//*                                 Moved GetEnergyLoss() and CalcQms()
//*                                 from TVKalDetector.
//*   2011/12/03  S.Aplin           Added new member: name 
//*                                 default value set to "TVMeasLayer"
//*                                 and corresponding member function
//*                                 TString GetName()
//*
//*************************************************************************

#include "TVMeasLayer.h"  // from KalTrackLib
#include "TKalTrack.h"    // from KalTrackLib
#include "TVTrack.h"      // from KalTrackLib

ClassImp(TVMeasLayer)

//_________________________________________________________________________
// ----------------------------------
//  Ctor
// ----------------------------------
TVMeasLayer::TVMeasLayer(TMaterial     &matIn,
                         TMaterial     &matOut,
                         Bool_t         isactive,
                         const Char_t    *name)
           : fMaterialInPtr(&matIn),
             fMaterialOutPtr(&matOut),
             fIndex(0),
             fIsActive(isactive),
             fname(name)
{
}

//_________________________________________________________________________
//  ----------------------------------
//   Utility Methods
//  ----------------------------------
//_________________________________________________________________________
// -----------------
//  GetEnergyLoss
// -----------------
//    returns energy loss.
//
Double_t TVMeasLayer::GetEnergyLoss(      Bool_t    isoutgoing,
                                    const TVTrack  &hel,
                                          Double_t  df) const
{
   Double_t cpa    = hel.GetKappa();
   Double_t tnl    = hel.GetTanLambda(); 
   Double_t tnl2   = tnl * tnl;
   Double_t tnl21  = 1. + tnl2;
   Double_t cslinv = TMath::Sqrt(tnl21);
   Double_t mom2   = tnl21 / (cpa * cpa);

   // -----------------------------------------
   // Bethe-Bloch eq. (Physical Review D P195.)
   // -----------------------------------------
   static const Double_t kK   = 0.307075e-3;     // [GeV*cm^2]
   static const Double_t kMe  = 0.510998902e-3;  // electron mass [GeV]
   static const Double_t kMpi = 0.13957018;      // pion mass [GeV]

   TKalTrack *ktp  = static_cast<TKalTrack *>(TVKalSystem::GetCurInstancePtr());
   Double_t   mass = ktp ? ktp->GetMass() : kMpi;

   const TMaterial &mat = GetMaterial(isoutgoing);
   Double_t dnsty = mat.GetDensity();		// density
   Double_t A     = mat.GetA();                 // atomic mass
   Double_t Z     = mat.GetZ();                 // atomic number
   //Double_t I    = Z * 1.e-8;			// mean excitation energy [GeV]
   //Double_t I    = (2.4 +Z) * 1.e-8;		// mean excitation energy [GeV]
   Double_t I    = (9.76 * Z + 58.8 * TMath::Power(Z, -0.19)) * 1.e-9;
   Double_t hwp  = 28.816 * TMath::Sqrt(dnsty * Z/A) * 1.e-9;
   Double_t bg2  = mom2 / (mass * mass);
   Double_t gm2  = 1. + bg2;
   Double_t meM  = kMe / mass;
   Double_t x    = log10(TMath::Sqrt(bg2));
   Double_t C0   = - (2. * log(I/hwp) + 1.);
   Double_t a    = -C0/27.;
   Double_t del;
   if (x >= 3.)            del = 4.606 * x + C0;
   else if (0.<=x && x<3.) del = 4.606 * x + C0 + a * TMath::Power(3.-x, 3.);
   else                    del = 0.;
   Double_t tmax = 2.*kMe*bg2 / (1. + meM*(2.*TMath::Sqrt(gm2) + meM)); 
   Double_t dedx = kK * Z/A * gm2/bg2 * (0.5*log(2.*kMe*bg2*tmax / (I*I))
                 - bg2/gm2 - del);

   Double_t path = hel.IsInB()
                 ? TMath::Abs(hel.GetRho()*df)*cslinv
                 : TMath::Abs(df)*cslinv;

   //fg: switched from using cm to mm in KalTest - material (density) and energy still in GeV and cm
   path /= 10. ; 

   Double_t edep = dedx * dnsty * path;


   Double_t cpaa = TMath::Sqrt(tnl21 / (mom2 + edep
                 * (edep + 2. * TMath::Sqrt(mom2 + mass * mass))));
   Double_t dcpa = TMath::Abs(cpa) - cpaa;

   static const Bool_t kForward  = kTRUE;
   static const Bool_t kBackward = kFALSE;
   Bool_t isfwd = ((cpa > 0 && df < 0) || (cpa <= 0 && df > 0)) ? kForward : kBackward;
   return isfwd ? (cpa > 0 ? dcpa : -dcpa) : (cpa > 0 ? -dcpa : dcpa);
}

//_________________________________________________________________________
// -----------------
//  CalQms
// -----------------
//    calculates process noise matrix for multiple scattering with
//    thin layer approximation.
//
void TVMeasLayer::CalcQms(      Bool_t       isoutgoing,
                          const TVTrack     &hel,
                                Double_t     df,
                                TKalMatrix  &Qms) const
{
   Double_t cpa    = hel.GetKappa();
   Double_t tnl    = hel.GetTanLambda(); 
   Double_t tnl2   = tnl * tnl;
   Double_t tnl21  = 1. + tnl2;
   Double_t cpatnl = cpa * tnl;
   Double_t cslinv = TMath::Sqrt(tnl21);
   Double_t mom    = TMath::Abs(1. / cpa) * cslinv;

   static const Double_t kMpi = 0.13957018; // pion mass [GeV]
   TKalTrack *ktp  = static_cast<TKalTrack *>(TVKalSystem::GetCurInstancePtr());
   Double_t   mass = ktp ? ktp->GetMass() : kMpi;
   Double_t   beta = mom / TMath::Sqrt(mom * mom + mass * mass);

   const TMaterial &mat = GetMaterial(isoutgoing);
   Double_t x0inv = 1. / mat.GetRadLength();  // radiation length inverse

   // *Calculate sigma_ms0 =============================================
   static const Double_t kMS1  = 0.0136;
   static const Double_t kMS12 = kMS1 * kMS1;
   static const Double_t kMS2  = 0.038;
   
   Double_t path = hel.IsInB()
                 ? TMath::Abs(hel.GetRho()*df)*cslinv
                 : TMath::Abs(df)*cslinv;

   //fg: switched from using cm to mm in KalTest - material (density) and energy still in GeV and cm
   path /= 10. ; 

   Double_t xl   = path * x0inv;
   // ------------------------------------------------------------------
   // Very Crude Treatment!!
   Double_t tmp = 1. + kMS2 * TMath::Log(TMath::Max(1.e-4, xl));
   tmp /= (mom * beta);
   Double_t sgms2 = kMS12 * xl * tmp * tmp;
   // ------------------------------------------------------------------

   Qms(1,1) = sgms2 * tnl21;
   Qms(2,2) = sgms2 * cpatnl * cpatnl;
   Qms(2,4) = sgms2 * cpatnl * tnl21;
   Qms(4,2) = sgms2 * cpatnl * tnl21;
   Qms(4,4) = sgms2 * tnl21  * tnl21;
}
