//*************************************************************************
//* ====================
//*  TVKalSystem Class
//* ====================
//*
//* (Description)
//*   This is the base class of Kalman filter.
//* (Requires)
//* 	TObject
//* (Provides)
//* 	class TVKalSystem
//* (Update Recored)
//*   2003/09/30  K.Fujii	Original version.
//*   2009/06/18  K.Fujii       Implement inverse Kalman filter
//*
//*************************************************************************

#include <iostream>
#include <cstdlib>
#include "TVKalSystem.h"
#include "TVKalState.h"

//_____________________________________________________________________
//  ------------------------------
//  Base Class for measurement vector used by Kalman filter
//  ------------------------------
//
ClassImp(TVKalSystem)

TVKalSystem *TVKalSystem::fgCurInstancePtr = 0;

TVKalSystem::TVKalSystem(Int_t n) 
            :TObjArray(n),
             fCurSitePtr(0),
             fChi2(0.)
{
   if (!fgCurInstancePtr) fgCurInstancePtr = this;
}

TVKalSystem::~TVKalSystem() 
{
   if (this == fgCurInstancePtr) fgCurInstancePtr = 0;
}

//-------------------------------------------------------
// AddAndFilter 
//-------------------------------------------------------

Bool_t TVKalSystem::AddAndFilter(TVKalSite &next)
{
   SetCurInstancePtr(this);

   //
   // Propagate current state to the next site
   //

   GetState(TVKalSite::kFiltered).Propagate(next);
   
   //
   // Calculate new pull and gain matrix
   //

   if (next.Filter()) {
      //
      // Add this to the system if accepted.
      //

      Add(&next);
      fChi2 += next.GetDeltaChi2();
      return kTRUE;
   } else {
      return kFALSE; 
   }
}

//-------------------------------------------------------
// GetNDF
//-------------------------------------------------------

Int_t TVKalSystem::GetNDF(Bool_t self)
{
   Int_t ndf    = 0;
   Int_t nsites = GetEntries();
   for (Int_t isite=1; isite<nsites; isite++) {
       TVKalSite &site = *static_cast<TVKalSite *>(At(isite));
       if (!site.IsLocked()) ndf += site.GetDimension();
   }
   if (self) ndf -= GetCurSite().GetCurState().GetDimension();
   return ndf;
}

//-------------------------------------------------------
// SmoothBackTo 
//-------------------------------------------------------

void TVKalSystem::SmoothBackTo(Int_t k)
{
   TIter previous(this,kIterBackward);
   TIter cur     (this,kIterBackward);

   TVKalSite  *prePtr;
   TVKalSite  *curPtr = static_cast<TVKalSite *>(cur());
   TVKalState &cura   = curPtr->GetState(TVKalSite::kFiltered);
   TVKalState &scura  = curPtr->GetState(TVKalSite::kSmoothed); 
   if (!&scura) {
      curPtr->Add(&curPtr->CreateState(cura, cura.GetCovMat(),
                                       TVKalSite::kSmoothed));
   }

   while ((curPtr = static_cast<TVKalSite *>(cur())) && 
          (prePtr = static_cast<TVKalSite *>(previous()))) {
      curPtr->Smooth(*prePtr);
      fCurSitePtr = curPtr;
      if (IndexOf(curPtr) == k) break;
   }
}

//-------------------------------------------------------
// SmoothAll 
//-------------------------------------------------------

void TVKalSystem::SmoothAll()
{
   SmoothBackTo(0);
}

//-------------------------------------------------------
// InvFilter
//-------------------------------------------------------

void TVKalSystem::InvFilter(Int_t k)
{
   using namespace std;
   //
   // Check if site k exists
   //
   TVKalSite  *curPtr = static_cast<TVKalSite *>(At(k));
   if (!curPtr) {
      cerr << "::::: ERROR in TVKalSystem::InvFilter(k=" << k << ")"  << endl
           << "  Site " << k << " nonexistent! Abort!"
           << endl;
      ::abort();
   }
   //
   // Check if site k already smoothed
   //
   if (!&curPtr->GetState(TVKalSite::kSmoothed)) SmoothBackTo(k);
   //
   // Inverse filter site k
   //
   fCurSitePtr = curPtr;
   curPtr->InvFilter();
}
