#include "GearTPCHit.h"

#include <iomanip>
// #include <streamlog/streamlog.h>

namespace kaldet{

//_________________________________________________________________________
//  --------------
//  Ctors and Dtor
//  --------------
//
GearTPCHit::GearTPCHit(Int_t m)
        : TVTrackHit(m),
          fXXPtr(0),
          fHitPtr(0),
	  fVDrift(0)
{
}

GearTPCHit::GearTPCHit(const TVMeasLayer &ms,
                         Double_t       *x,
                         Double_t       *dx,
                   const TVector3       &xx,
                         Double_t        b,
		         Double_t        v,
                         Int_t           m)
        : TVTrackHit(ms, x, dx, b, m),
	  fXXPtr(&xx),
          fHitPtr(0),
	  fVDrift(v)
{
}

GearTPCHit::GearTPCHit(const TVMeasLayer &ms,
		       Double_t       *x,
		       Double_t       *dx,
		       const void           *hitp,
		       Double_t        b,
		       Double_t        v,
		       Int_t           m)
        : TVTrackHit(ms, x, dx, b, m),
          fXXPtr(0),
          fHitPtr(hitp),
	  fVDrift(v)
{
}

GearTPCHit::~GearTPCHit()
{
}

//_________________________________________________________________________
//  ----------------
//  Compare two Hits
//  ----------------
//
Int_t GearTPCHit::Compare(const TObject *obj) const
{
  Double_t me  = GetSortingPolicy();
  Double_t you = dynamic_cast<const GearTPCHit *>(obj)->GetSortingPolicy();
  if (you == 0) 
  {
    // streamlog_out(ERROR) << "Cannot compare GearTPCHit to something which is not a"
    //     		 << " GearTPCHit" << std::endl;
    throw std::bad_cast();
  }
  return me < you ? -1 : (me > you ? +1 : 0);
}

Double_t GearTPCHit::GetSortingPolicy() const
{
  return GetMeasLayer().HitToXv(*this).Mag();
}

Bool_t GearTPCHit::IsSortable() const
{
  return kTRUE;
}

}//namespace kaldet
