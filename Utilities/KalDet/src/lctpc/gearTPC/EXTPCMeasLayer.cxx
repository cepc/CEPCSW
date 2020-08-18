#include "EXTPCMeasLayer.h"
#include "EXTPCKalDetector.h"
#include "EXTPCHit.h"

#include <TRandom.h>
// #include <streamlog/streamlog.h>

//ClassImp(EXTPCMeasLayer)

EXTPCMeasLayer::EXTPCMeasLayer(TMaterial &min,
			       TMaterial &mout,
			       Int_t      module,
			       Int_t      row,
			       Double_t   r0,
			       Double_t   lhalf,
			       TVector3   xc,
			       Bool_t     isPerfect,
			       Bool_t     isActive,
			       Double_t   sigmaX0,
			       Double_t   sigmaX1,
			       Double_t   sigmaZ0,
			       Double_t   sigmaZ1)
  :  kaldet::GearTPCMeasLayer(min, mout, module, row, isPerfect, isActive,
			      sigmaX0, sigmaX1, sigmaZ0, sigmaZ1),
     TCylinder(r0, lhalf, xc.X(), xc.Y(), xc.Z())

{
}

EXTPCMeasLayer::~EXTPCMeasLayer(){}

Int_t EXTPCMeasLayer::GetModuleID() const
{
  return fModuleRows.begin()->first;
}

Int_t EXTPCMeasLayer::GetLayerID() const
{
  return fModuleRows.begin()->second;
}

TKalMatrix EXTPCMeasLayer::XvToMv(const TVector3 &xv, Int_t side) const
{
  return XvToMv(xv); 
}

void EXTPCMeasLayer::ProcessHit(const TVector3  &xx,
				  TObjArray &hits) const
{
  TKalMatrix h    = XvToMv(xx);
  Double_t   rphi = h(0, 0);
  Double_t   d    = h(1, 0);

  Double_t dx = GetSigmaX(d);
  Double_t dz = GetSigmaZ(d);
  rphi += gRandom->Gaus(0., dx);  // smearing rphi
  d    += gRandom->Gaus(0., dz);  // smearing drift distance

  Double_t meas [2];
  Double_t dmeas[2];
  meas [0] = rphi;
  meas [1] = d;
  dmeas[0] = dx;
  dmeas[1] = dz;

  Double_t b;
  Double_t vDrift;
  try{
    EXTPCKalDetector const & detector = dynamic_cast<EXTPCKalDetector const &>(GetParent(false));
    b = detector.GetBfield();
    vDrift = detector.GetVdrift();
  }
  catch( std::bad_cast & )
  {
    // streamlog_out(ERROR) << "EXTPCMeasLayer::ProcessHit: Error: Parent is not an EXTPCKalDetector."
    //     		 << std::endl;
    // streamlog_out(ERROR) << "This function is for backward compatibility only."
    //     		 <<" Please do not use it when directly instantiating a GearTPCKalDetector."
    //     		 << std::endl;    
    throw;
  }
  
  hits.Add(new EXTPCHit(*this, meas, dmeas, 1 /*side, not functional anyway*/,
			vDrift, xx, b ));
}
