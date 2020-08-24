#include "EXTPCHit.h"

EXTPCHit::EXTPCHit(Int_t m) : kaldet::GearTPCCylinderHit(m) , fSide(0) {}

EXTPCHit::EXTPCHit(const TVMeasLayer &ms,
		   Double_t       *x,
		   Double_t       *dx,
		   Int_t           side,
		   Double_t        v,
		   const TVector3       &xx,
		   Double_t        b,
		   Int_t           m)
  : kaldet::GearTPCCylinderHit(ms, x, dx, xx, b, v, m), fSide(side) {}

EXTPCHit::EXTPCHit(const TVMeasLayer &ms,
		   Double_t       *x,
		   Double_t       *dx,
		   Int_t           side,
		   Double_t        v,
		   const void           *hitp,
		   Double_t        b,
		   Int_t           m)
  : kaldet::GearTPCCylinderHit(ms, x, dx, hitp, b, v, m), fSide(side) {}

EXTPCHit::~EXTPCHit()
{}

