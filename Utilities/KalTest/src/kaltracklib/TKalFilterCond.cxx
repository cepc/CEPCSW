//*************************************************************************
//* =====================
//*  TKalFilterCond Class
//* =====================
//*
//* (Description)
//*   A class to specify filter conditions used in Kalman filter.
//* (Requires)
//* (Provides)
//* 	class TKalFilterCond
//* (Update Recored)
//*   2010/04/06  K.Fujii        Original Version.
//*
//*************************************************************************

#include "TKalFilterCond.h"
#include "TKalTrackSite.h"
#include <iostream>

//_____________________________________________________________________
//  ------------------------------
//   Filter condition class
//  ------------------------------

ClassImp(TKalFilterCond)

Bool_t TKalFilterCond::IsAccepted(const TKalTrackSite &site)
{
#if 0
   // return kTRUE if this site is acceptable
   using namespace std;
   Double_t delchi2 = site.GetDeltaChi2();
   if (delchi2 > 25.) {
      cerr << ">>>> TKalFilterCond::IsAccepted >>>>>>>>>>>>> " << endl
           << " Too big chi2 increment!!! " << endl;
      site.DebugPrint();
      return kFALSE;
   }
#else
   return kTRUE;
#endif
}
