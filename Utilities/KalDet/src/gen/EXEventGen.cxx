#include "EXEventGen.h"
//#include "EXVKalDetector.h"
//#include <EXVMeasLayer.h>
#include "EXTPCKalDetector.h"
#include <EXTPCMeasLayer.h>
#include <kaltest/TPlane.h>
#include <TRandom.h>
#include <TMath.h>

#include <iostream>
#include <typeinfo> // needed for std::bad_cast exception

//-----------------------------------
// Track Parameters
//-----------------------------------

#define __DR__     0.
#define __DZ__     0.

Double_t EXEventGen::fgT0 = 0.; // [nsec]

THelicalTrack EXEventGen::GenerateHelix(Double_t pt,
                                        Double_t cosmin,
                                        Double_t cosmax,
					Double_t phimin,
					Double_t phimax,
					TVector3 xv0)
{
   // ---------------------------
   //  Generate a helical track
   // ---------------------------

   Double_t dr  = __DR__;
   Double_t fi0 = gRandom->Uniform(phimin,phimax);
   Double_t cpa = 1. / pt;
   Double_t dz  = __DZ__;
   Double_t cs  = gRandom->Uniform(cosmin, cosmax);
   Double_t tnl = cs / TMath::Sqrt((1-cs)*(1+cs)); 
   Double_t x0  = xv0.X();
   Double_t y0  = xv0.Y();
   Double_t z0  = xv0.Z();

   EXTPCMeasLayer * measLayer = dynamic_cast<EXTPCMeasLayer *>(fCradlePtr->At(0));
   if (measLayer==0)
   {
     std::cerr << "EXEventGen::GenerateHelix: cast of fCradlePtr->At(0) to EXVMeasLayer * failed"
	       << std::endl;
     std::cerr << "Derive your measurement layer from EXTPCMeasLayer instead of "
	       << "TVMeasLayer to use EXEventGen." <<  std::endl;
   
     throw std::bad_cast();
   }

   EXTPCKalDetector const & kalDetector = dynamic_cast<const EXTPCKalDetector &>(measLayer->GetParent(kFALSE));
   if (measLayer==0)
   {
     std::cerr << "EXEventGen::GenerateHelix: cast of (measLayer->GetParent(kFALSE) to EXTPCKalDetector * failed" << std::endl;
     std::cerr << "Derive your detector from EXTPCKalDetector instead of "
	       << "TVKalDetector to use EXEventGen." <<  std::endl;
     throw std::bad_cast();
   }

   Double_t b   = kalDetector.GetBfield();

   return THelicalTrack(dr,fi0,cpa,dz,tnl,x0,y0,z0,b);
}

void EXEventGen::Swim(THelicalTrack &heltrk)
{
   // ---------------------------
   //  Swim track and Make hits
   // ---------------------------

   // killenb: I have no idea what dfi is (delta phi?), but dividing the sorting policy by a 
   // radius does not make sense. The sorting polily used to return  GetR() + GetXc().X(),
   // so let's use this. It will onyl work for the LP1, but this class is deprecated anyway.
   TCylinder * firstMeasCylinder = dynamic_cast<TCylinder *>(fCradlePtr->At(0));
   if (firstMeasCylinder==0)
   {
     std::cerr << "Cannot cast first object in the cradle to TCylinder*"<< std::endl;
     std::cerr << "EXEventGen only works for cylindrical measuremnt layers, sorry."<< std::endl;
     throw std::bad_cast();
   }
   Double_t dfi       = - (firstMeasCylinder->GetR() + firstMeasCylinder->GetXc().X())
                         / heltrk.GetRho();

   Int_t    nLayers   = fCradlePtr->GetEntries();
   Int_t    dLayer      = 1;
   Double_t dfisum    = 0.;

   for (Int_t layer = 0; layer >= 0; layer += dLayer) { // loop over layers
      // change direction if it starts looping back
      if (layer == nLayers - 1) dLayer = -1;

      // FIXME: This assumes that the measurement layer implementation is derived from
      // TVSurface. This is not ensured by the interface!
      EXTPCMeasLayer const * measurementLayer    
	= dynamic_cast<EXTPCMeasLayer const *>(fCradlePtr->At(layer));
      //check that the object is a measurement layer. It better should be...
      if (measurementLayer==0)
      {
	std::cerr << "EXEventGen::Swim(): fCradlePtr->At(layer) is not an EXVMeasLayer" << std::endl;
	throw std::bad_cast();	
      }
      
      TVSurface const * measurementSurface  = dynamic_cast<TVSurface   const *>(fCradlePtr->At(layer));
      //check that the object is a surface. This part is not ensured by the interface.
      if (measurementSurface==0)
      {
	std::cerr << "EXEventGen::Swim(): fCradlePtr->At(layer) is not a TVSurface" << std::endl;
	throw std::bad_cast();	
      }


      TVector3 xx;
      Double_t dfis = dfi;
      if (!measurementSurface->CalcXingPointWith(heltrk,xx,dfi,1)
       || TMath::Abs(dfi) > TMath::Pi()
       || TMath::Abs(dfi + dfisum) > TMath::TwoPi()) {
         dfi = dfis;
         continue;
      }
      // should use the material behind the surface since dfi is measured 
      // from the last point to the current surface
      // Bool_t   dir    = dKayer < 0 ? kTRUE : kFALSE;

      dfisum += dfi;

      heltrk.MoveTo(xx,dfi);	// move pivot to current hit

      // killenb: seriously: check if the detector is powered?
      if (measurementLayer->IsActive() ) 
	//  This has been kickes as true was hard coded anyway
	// && dynamic_cast<const EXVKalDetector &>(measurementLayer->GetParent(kFALSE)).IsPowerOn()) 
      {
	measurementLayer->ProcessHit(xx, *fHitBufPtr); // create hit point
      }
      if (layer == nLayers - 1) break;
   }
}
