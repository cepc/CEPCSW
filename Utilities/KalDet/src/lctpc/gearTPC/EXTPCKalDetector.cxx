#include "EXTPCKalDetector.h"
#include <kaltest/TKalDetCradle.h>
#include <gear/GearMgr.h>
#include <gear/BField.h>
// #include <streamlog/streamlog.h>
#include <stdlib.h>
// global constants from Marlin, used for the global pointer to the GearMgr
// KILLENB: Worst style I have ever seen. Unfortunately needed to implement
// the static (stupid!!) GetInstance function. I will kick this out as soon
// as I have the comparison with the old code.
//#include <marlin/Global.h>


Double_t EXTPCKalDetector::fgVdrift = 76.e-3; // [mm/nsec]
EXTPCKalDetector * EXTPCKalDetector::fgInstance=0;

//ClassImp(EXTPCKalDetector)

EXTPCKalDetector::EXTPCKalDetector(const gear::GearMgr& gearMgr)
  : kaldet::GearTPCKalDetector(gearMgr)
{
  gear::BField const & bField = gearMgr.getBField();
  // get the BField at 0,0,0. Check that there are no transverse components
  // FIXME: Event though there are no transverse components at 0, 0, 0 does not mean
  // there are no transverse components somewhere else.
  gear::Vector3D bFieldVector =  bField.at(gear::Vector3D(0., 0., 0.));
  if (bFieldVector[0]!=0 || bFieldVector[1]!=0)
  {
    // streamlog_out(ERROR) << "B-field has transverse components."
    //     		 << " EXTPCKalDetector only works with homogeneous B-field"
    //     		 << " in z direction" << std::endl;
    throw gear::Exception("Magnetic field not in z direction.");
  }
  
  fBField = bFieldVector[2];
}

EXTPCKalDetector::~EXTPCKalDetector()
{}

EXTPCKalDetector * EXTPCKalDetector::GetInstance()
{
  //gear::GearMgr const * gearMgr = marlin::Global::GEAR;

  if (!fgInstance) {
    //fgInstance = new EXTPCKalDetector(*gearMgr);
    std::cout << "need to fix, no gearMgr...." << std::endl;
    exit(1);
    // Let's leak some memory. This cradle can never be deleted because we loose the pointer.
    // But it cannot be deleted anyway because the Layers in the KalDetector have the cradle as
    // parent, so deleting the cradle voids the detector.
    // But you also cannot delete the KalDetector, because there are layers pointed to by the cradle.
    // Stupid design.
    TKalDetCradle * lp1tpc = new TKalDetCradle;
    lp1tpc->Install(*fgInstance);
    lp1tpc->Close();
    lp1tpc->Sort();

     //std::cout << "# of layers: " << lp1tpc.GetEntries()  << std::endl;
  }
  return fgInstance;
  
}

 Double_t EXTPCKalDetector::GetBfield()
 {
   if (fgInstance)
     return fgInstance->fBField;
   else
   {
     // streamlog_out(ERROR) << " EXTPCKalDetector::GetBField() called without instantiating the object"
     //    		  << std::endl;
     throw std::exception();
     return 0;
   }
 }
