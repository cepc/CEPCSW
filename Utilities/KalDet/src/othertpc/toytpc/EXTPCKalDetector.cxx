//*************************************************************************
//* ========================
//*  EXTPCKalDetector Class
//* ========================
//*
//* (Description)
//*   User defined detector class
//* (Requires)
//*     EXVKalDetector
//* (Provides)
//*     class EXTPCKalDetector
//* (Update Recored)
//*   2009/11/23  K.Ikematsu   Derived from KalTest/examples/kaltest/
//*                                         hybrid/tpc/EXTPCKalDetector.cxx
//*
//* $Id: EXTPCKalDetector.cxx,v 1.1 2010-03-11 15:07:01 fujiik Exp $
//*************************************************************************
//
// STL
#include <vector>

// GEAR
#include "gear/GEAR.h"
#include "gear/TPCParameters.h"
#include "gear/PadRowLayout2D.h"
#include "gearimpl/TPCModuleImpl.h"
#include "gearxml/GearXML.h"

// global constants from Marlin, used for the global pointer to the GearMgr
//#include <marlin/Global.h>

#include "EXTPCKalDetector.h"
#include "EXTPCMeasLayer.h"
#include "EXTPCHit.h"
#include "kaltest/TKalDetCradle.h"

// ROOT
#include "TRandom.h"
#include "TMath.h"
#include "TTUBE.h"
#include "TNode.h"
#include "TVirtualPad.h"

EXTPCKalDetector * EXTPCKalDetector::fgInstance = 0;
Double_t EXTPCKalDetector::fgVdrift = 7.6e-3;


ClassImp(EXTPCKalDetector)

// definition of gas parameters for mixture of Ar + CF4 + C4H10
Double_t Ar = 0.95;  // fraction of Argon in chamber gas
Double_t CF = 0.043; // fraction of CF4 in chamber gas
Double_t CH = 0.007; // fraction of C4H10 in chamber gas

EXTPCKalDetector::EXTPCKalDetector(Int_t m)
                : EXVKalDetector(m),
                  fNodePtr(0)
{
  Double_t A, Z, density, radlen;

  A       = 39.948 * Ar + 17.601 * CF + 4.152 * CH;
  Z       = 18 * Ar + 8.4 * CF + 2.4 * CH;
  density = (1.784 * Ar + 3.692 * CF + 2.65 * CH)/1000;
  radlen  = 1.196e4*2; // has to be checked, copied from 0.9 Ar, 0.1 C4H10 lines
  TMaterial &gas = *new TMaterial("TPCGas", "", A, Z, density, radlen, 0.);
#if 1
   static const Int_t    nlayers   = 24;            // number of layer
   static const Double_t kmm2cm    = 0.1;
   static const Double_t lhalf     = 600. * kmm2cm; //length
   // diffusion coefficients transversal (x) and longitudinal (z) to magnetic field
   static const Double_t neff      = 22.7;
#if 0
   static const Double_t sigmax0   = 38.3e-4;
   static const Double_t sigmax1   = 101.5e-4 /TMath::Sqrt(10.)/TMath::Sqrt(neff);
   static const Double_t sigmaz0   = 500.e-4;
   static const Double_t sigmaz1   = 154.e-4/TMath::Sqrt(10.)/TMath::Sqrt(neff);
#else
   static const Double_t sigmax0   = 50e-4;
   static const Double_t sigmax1   = 0.;
   static const Double_t sigmaz0   = 500.e-4;
   static const Double_t sigmaz1   = 0.;
#endif
                         fgVdrift  = 50.e-3 * kmm2cm;
#else

  gear::TPCParameters const &theTPCParameters
    = marlin::Global::GEAR->getTPCParameters();

  Int_t nmodules = theTPCParameters.getNModules();

  std::vector<const gear::TPCModule *> modules;

  for (Int_t i = 0; i < nmodules; i++) {
    modules.push_back(&theTPCParameters.getModule(i));
  }

  static const Double_t kmm2cm = 0.1;

  static const Double_t lhalf
    = theTPCParameters.getMaxDriftLength() * kmm2cm;     // half length
  static const Int_t    nrows = modules[0]->getNRows();  // # of pad rows
  ///// FIXME: temporary treatment /////////////////////////
  static const Int_t    nlayers = nrows * 3;  // # of layers
  //////////////////////////////////////////////////////////
  static const Double_t neff    = 22.7;
  static const Double_t sigmax0 = 38.3e-4;
  static const Double_t sigmax1 = 101.5e-4 / TMath::Sqrt(10.) / TMath::Sqrt(neff);
  static const Double_t sigmaz0 = 500.e-4;
  static const Double_t sigmaz1 = 154.e-4  / TMath::Sqrt(10.) / TMath::Sqrt(neff);
#endif
  Bool_t active = EXTPCMeasLayer::kActive;

  // create measurement layers
  Int_t module = 1;
  for (Int_t layer = 0; layer < nlayers; layer++) {
     Double_t XMin  = -32. * kmm2cm;
     Double_t XMax  = 32. * kmm2cm;
     Double_t y0    = (4.*layer-46.) * kmm2cm;
     Add(new EXTPCMeasLayer(gas, gas, y0, XMin, XMax, lhalf, sigmax0, sigmax1, sigmaz0, sigmaz1, active, layer, module));
  }
  SetOwner();
}

EXTPCKalDetector::~EXTPCKalDetector()
{
}

EXTPCKalDetector * EXTPCKalDetector::GetInstance()
{
  if (!fgInstance) {
     fgInstance = new EXTPCKalDetector;
     TKalDetCradle & toydet = * new TKalDetCradle;
     toydet.Install(*fgInstance);
     toydet.Close();
     toydet.Sort();
  }
  return fgInstance;
}

void EXTPCKalDetector::Draw(Int_t color, const Char_t *opt)
{
  if (! gPad) return;
  TNode *nodep = GetNodePtr();
  nodep->cd();
  //EXVKalDetector::Draw(color, opt);
  TAttDrawable::Draw(color, opt);
}
