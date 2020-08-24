#include "GearTPCMeasLayer.h"

#include <TMath.h>

#include <iostream>
// #include <streamlog/streamlog.h>

#include <gear/GEAR.h>

namespace kaldet
{


GearTPCMeasLayer::GearTPCMeasLayer(TMaterial &min,
				   TMaterial &mout,
				   Int_t      module,
				   Int_t      row,
				   Bool_t     isPerfect,
				   Bool_t     isActive,
				   Double_t   sigmax0,
				   Double_t   sigmax1,
				   Double_t   sigmaz0,
				   Double_t   sigmaz1)
  : TVMeasLayer(min, mout, isActive),
    fSigmaX0(sigmax0),
    fSigmaX1(sigmax1),
    fSigmaZ0(sigmaz0),
    fSigmaZ1(sigmaz1),
    fIsPerfect(isPerfect)
{
  fModuleRows.insert( std::pair<int, int>(module, row) );
}

GearTPCMeasLayer::~GearTPCMeasLayer()
{
}

Bool_t GearTPCMeasLayer::IsPerfect() const
{
  return fIsPerfect;
}

std::set< std::pair <int, int> > const & GearTPCMeasLayer::GetModuleRows() const
{
  return fModuleRows;
}


void GearTPCMeasLayer::AddModuleRow(int module, int row)
{
  if (fIsPerfect)
  {
    fModuleRows.insert( std::pair<int, int>(module, row) );
  }
  else
  {
    // streamlog_out(ERROR) << "You can only add additional modules to perfect layers." <<std::endl;
    throw gear::Exception("GearTPCMeasLayer is not declared perfect.");
  }
}

Double_t GearTPCMeasLayer::GetSigmaX(Double_t zdrift) const
{
  return TMath::Sqrt(fSigmaX0 * fSigmaX0 + fSigmaX1 * fSigmaX1 * zdrift);
}

Double_t GearTPCMeasLayer::GetSigmaZ(Double_t zdrift) const
{
  return TMath::Sqrt(fSigmaZ0 * fSigmaZ0 + fSigmaZ1 * fSigmaZ1 * zdrift);
}

}//namespace kaldet
