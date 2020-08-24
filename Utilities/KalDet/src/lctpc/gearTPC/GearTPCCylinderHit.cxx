#include "GearTPCCylinderHit.h"
#include "GearTPCCylinderMeasLayer.h"
#include <TMath.h>

#include <iostream>
#include <iomanip>
// #include <streamlog/streamlog.h>

using namespace std;

namespace kaldet{

//_________________________________________________________________________
//  --------------
//  Ctors and Dtor
//  --------------
//
GearTPCCylinderHit::GearTPCCylinderHit(Int_t m)
        : GearTPCHit(m)
{
}

GearTPCCylinderHit::GearTPCCylinderHit(const TVMeasLayer &ms,
			      Double_t       *x,
			      Double_t       *dx,
			      const TVector3       &xx,
			      Double_t        b,
			      Double_t        v,
			      Int_t           m)
  : GearTPCHit(ms, x, dx, xx, b, v, m)
{
}

GearTPCCylinderHit::GearTPCCylinderHit(const TVMeasLayer &ms,
			       Double_t       *x,
			       Double_t       *dx,
			       const void           *hitp,
			       Double_t        b,
			       Double_t        v,
			       Int_t           m)
  : GearTPCHit(ms, x, dx, hitp, b, v, m)
{
}

GearTPCCylinderHit::~GearTPCCylinderHit()
{
}

//_________________________________________________________________________
//  --------------------------------
//  Implementation of public methods
//  --------------------------------
//
TKalMatrix GearTPCCylinderHit::XvToMv(const TVector3 &xv, Double_t t0) const
{
  // KILLENB What does this do? The ild version just returns measLayer->XvToMv()

  TKalMatrix h(GetDimension());
  Double_t   r;
  try{
    const GearTPCCylinderMeasLayer &ms
      = dynamic_cast<const GearTPCCylinderMeasLayer &>(GetMeasLayer());
    h = ms.XvToMv(xv);
    r = ms.GetR();
  }
  catch(std::bad_cast &)
  {
    // streamlog_out(ERROR) << "GearTPCHit::XvToMv: Measurement layer is not a GearTPCCylinderMeasLayer\n"
    //     		 << "  The current implementation only works for cylindrical layer, sorry."
    //     		 << std::endl;
    throw;
  }


  Double_t   phih = h(0, 0) / r;
  Double_t   phim = (*this)(0, 0) / r;
  Double_t   dphi = phih - phim;

  static const Double_t kPi    = TMath::Pi();
  static const Double_t kTwoPi = 2 * kPi;

  while (dphi < -kPi) dphi += kTwoPi;
  while (dphi >  kPi) dphi -= kTwoPi;

  h(0, 0)  = r * (phim + dphi);
  h(1, 0) += fVDrift * t0;

  return h;
}

void GearTPCCylinderHit::DebugPrint(Option_t *) const
{
  cerr << "------------------- Site Info -------------------------" << endl;

  for (Int_t i = 0; i < GetDimension(); i++) {
    Double_t x  = (*this)(i, 0);
    Double_t dx = (*this)(i, 1);
    cerr << " x[" << i << "] = " << setw(8) << setprecision(5) << x
         << "    "
         << "dx[" << i << "] = " << setw(6) << setprecision(2) << dx
         << setprecision(7)
         << resetiosflags(ios::showpoint)
         << endl;
  }
  cerr << " r    = " << setw(8)
       << dynamic_cast<const GearTPCCylinderMeasLayer &>(GetMeasLayer()).GetR() << endl;
  cerr << "-------------------------------------------------------"  << endl;
}


}//namespace kaldet
