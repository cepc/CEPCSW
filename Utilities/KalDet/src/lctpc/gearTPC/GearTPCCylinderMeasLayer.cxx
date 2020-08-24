#include "GearTPCCylinderMeasLayer.h"
#include "GearTPCHit.h"
#include "GearTPCKalDetector.h"
#include "GearTPCCylinderHit.h"

#include <TMath.h>

#include <iostream>
// #include <streamlog/streamlog.h>

#include <gear/GEAR.h>

namespace kaldet
{


GearTPCCylinderMeasLayer::GearTPCCylinderMeasLayer(TMaterial &min,
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
						   Double_t   sigmaZ1,
						   Double_t   phimin,
						   Double_t   phimax)
  : EXTPCMeasLayer(min, mout, module, row, r0, lhalf, xc ,isPerfect, 
		   isActive, sigmaX0, sigmaX1, sigmaZ0, sigmaZ1),
    /* this is the original code which should be reactivated once the EXTPCMeasLayer is phased out:
       : GearTPCMeasLayer(min, mout, module, isPerfect, isActive,
                          sigmaX0, sigmaX1, sigmaZ0, sigmaZ1),
         TCylinder(r0, lhalf, xc.X(), xc.Y(), xc.Z()),
    */
                fPhiMin(phimin),
                fPhiMax(phimax)
{
  //FIXME: As the handling of cylinder segments is not defined yet we force the layer to be
  //a perfect cylinder
  fIsPerfect = true;

  // for a full cylinder phi min and max do not make sense. Give a warning if they are used.
  if (fIsPerfect)
  {
    if (phimin!=-TMath::Pi())
    {
      // streamlog_out(WARNING) << "GearTPCCylinderMeasLayer: Ignoring input parameter phimin."
      //   		     << " The current implementation is a full cylinder." << std::endl;
      phimin=-TMath::Pi();
    }
    if (phimax!=TMath::Pi())
    {
      // streamlog_out(WARNING) << "GearTPCCylinderMeasLayer: Ignoring input parameter phimax."
      //   		     << " The current implementation is a full cylinder." << std::endl;
      phimax=TMath::Pi();
    }
  }
}

GearTPCCylinderMeasLayer::~GearTPCCylinderMeasLayer()
{
}

TKalMatrix GearTPCCylinderMeasLayer::XvToMv(const TVector3 &xv) const
{
  TVector3 xxv = xv - GetXc();

  Double_t phi = TMath::ATan2(xxv.Y(), xxv.X());// - fPhiMin;

  static Double_t kPi    = TMath::Pi();
  static Double_t kTwoPi = 2 * kPi;
  while (phi < -kPi) phi += kTwoPi;
  while (phi >  kPi) phi -= kTwoPi;

  // Calculate hit coordinate information:
  //   mv(0, 0) = r * phi
  //     (1, 0) = drift distance

  TKalMatrix mv(kMdim, 1);
  mv(0, 0) = GetR() * phi;

  mv(1, 0) = xxv.Z();

  return mv;
}

TKalMatrix GearTPCCylinderMeasLayer::XvToMv(const TVTrackHit &vhit,
                                  const TVector3   &xv) const
{
  return XvToMv(xv);
}

TVector3 GearTPCCylinderMeasLayer::HitToXv(const TVTrackHit &vhit) const
{
  //  const EXTPCHit &hit = dynamic_cast<const EXTPCHit &>(vhit);


  Double_t phi = vhit(0, 0) / GetR();// + fPhiMin;

  Double_t z   = vhit(1, 0);

  Double_t x   = GetR() * TMath::Cos(phi) + GetXc().X();
  Double_t y   = GetR() * TMath::Sin(phi) + GetXc().Y();

  return TVector3(x, y, z);
}

void GearTPCCylinderMeasLayer::CalcDhDa(const TVTrackHit & vhit,
                              const TVector3   &xxv,
                              const TKalMatrix &dxphiada,
                                    TKalMatrix &H) const
{
  double vDrift;
  try{
    const GearTPCHit &hit = dynamic_cast<const GearTPCHit &>(vhit);
    vDrift = hit.GetVdrift();
  }
  catch(std::bad_cast &)
  {
    // streamlog_out(ERROR) << "GearTPCCylinderMeasLayer::CalcDhDa :"
    //     		 << "Cannot cast incoming TVTrackHit to GearTPCHit!"<< std::endl;
    throw;
  }

  // Calculate
  //    H = (@h/@a) = (@phi/@a, @z/@a)^t
  //  where
  //        h(a) = (phi, z)^t: expected meas vector
  //        a = (drho, phi0, kappa, dz, tanl, t0)
  //

  Int_t sdim = H.GetNcols();
  Int_t hdim = TMath::Max(5, sdim - 1);

  TVector3 xxvc = xxv - GetXc();
  Double_t xv   = xxvc.X();
  Double_t yv   = xxvc.Y();
  Double_t xxyy = xv * xv + yv * yv;

  // Set H = (@h/@a) = (@d/@a, @z/@a)^t

  for (Int_t i = 0; i < hdim; i++) {
    H(0, i)  = - (yv / xxyy) * dxphiada(0, i)
               + (xv / xxyy) * dxphiada(1, i);
    H(0, i) *= GetR();

    H(1, i)  = dxphiada(2, i);
  }

  if (sdim == 6) {
    H(0, sdim - 1) = 0.;

    // KILLENB I don't understand what this does, but I removed vDrift, so I set this to 0.
    H(1, sdim - 1) = - vDrift;
  }
}

Double_t GearTPCCylinderMeasLayer::GetSortingPolicy() const
{
   // The sorting policy (copied from the header file):
  // The layers are first sorted by radius + offset. This offset is only
  // useful for segments of a cylinder, like the LP1.
  // As offsets in this case can be positive or negative, but only make sense in one 
  // direction (you need a continuous number), we only allow offsets in X.
  // This should not be too much of a problem, you should be able to rotate your coordinates
  // so the offset is in X. If not you have to extend the sorting policy. (Please thake
  // care not to reduce versatility when doing so. You might want to implement your own class?)
  // 
  // For equal radii  + offset the layers are sorted by moduleID. As we have to squeeze this 
  // information into only one number, we multiply the radius + offset by 1e9 and add the moduleID.
  // A double has a precision of 53 bits, which is 15.9 digits. So the radius can be up to 1e6.9 mm
  // without causing the last digit of the the ModuleID to be cut, and for up to 1000 modules the
  // layers can be distinguished down to 1 nm without the two numbers mixing, or down to 1 micron
  // with up to 1.000.000 modules.
  
  // The additional sorting by module is intended for cylinder segments. Here only one module/row
  // per layer is allowed, so we just take the first entry in the set. In case of a perfect layer
  // it does not matter because there should only be one layer at this radius, so the sort order
  // should not be affected by adding an arbitrary module ID (as long as the module ID is < 1e6, as 
  // described above).
  
  // check that the offset is onyl in X
  if ( (GetXc().Y()!=0) || (GetXc().Z()!=0) )
  {
    throw gear::NotImplementedException("Offset is only allowed in X in the current implementation");
  }

  int moduleID = fModuleRows.begin()->first;
  // give a warning in case of very large module IDs. The sorting policy might have to be adapted
  if( moduleID > 1e6 )
  {
    // streamlog_out(WARNING) << "GearTPCCylinderMeasLayer::GetSortingPolicy(): "
    //     		   << "Very large module ID " << moduleID << " found. "
    //     		   << "This might compromise the sorting policy."
    //     		   << std::endl;
  }

  // FIXME: Do we have to check for a max r-offsetX?
  return (GetR() + GetXc().X())*1e9 + moduleID;
}

GearTPCHit * GearTPCCylinderMeasLayer::createHit(Double_t * meas,
						 Double_t * dmeas,
						 void * hitPointer, 
						 Double_t bField,
						 Double_t vDrift,
						 Int_t           m) const
{
  return new GearTPCCylinderHit(*this, meas, dmeas, hitPointer,
				bField, vDrift, m);
}


}//namespace kaldet
