#ifndef TVKALSYSTEM_H
#define TVKALSYSTEM_H
//*************************************************************************
//* ===================
//*  TVKalSystem Class
//* ===================
//*
//* (Description)
//*   Base class for Kalman filtering class
//* (Requires)
//* 	TObjArray
//* (Provides)
//* 	class TVKalSystem
//* (Update Recored)
//*   2003/09/30  K.Fujii	Original version.
//*   2005/08/25  A.Yamaguchi	Added fgCurInstancePtr and its getter & setter.
//*   2009/06/18  K.Fujii       Implement inverse Kalman filter
//*
//*************************************************************************

#include "TObjArray.h"
#include "TVKalSite.h"

//_____________________________________________________________________
//  ------------------------------
//  Kalman Filtering class
//  ------------------------------
//
class TKalMatrix;

class TVKalSystem : public TObjArray {
friend class TVKalSite;
public:

   // Ctors and Dtor

   TVKalSystem(Int_t n = 1);
   virtual ~TVKalSystem();

   // Utility methods

   virtual Bool_t AddAndFilter(TVKalSite &next);
   virtual void   SmoothBackTo(Int_t k);
   virtual void   SmoothAll();
   virtual void   InvFilter(Int_t k);

   inline  void   Add(TObject *obj);

   // Getters

   inline virtual TVKalSite  & GetCurSite()     { return *fCurSitePtr; }
   inline virtual TVKalState & GetState(TVKalSite::EStType t) 
                                   { return fCurSitePtr->GetState(t); }
   inline virtual Double_t     GetChi2() { return fChi2; }
          virtual Int_t        GetNDF (Bool_t self = kTRUE);
   
   static         TVKalSystem *GetCurInstancePtr() { return fgCurInstancePtr; }

   // Setters

private:
   static void SetCurInstancePtr(TVKalSystem *ksp) { fgCurInstancePtr = ksp; }

private:
   TVKalSite   *fCurSitePtr;  // pointer to current site
   Double_t     fChi2;        // current total chi2

   static TVKalSystem *fgCurInstancePtr;  //! currently active instance
   
   ClassDef(TVKalSystem,1)  // Base class for Kalman Filter
};

//=======================================================
// inline functions
//=======================================================

void TVKalSystem::Add(TObject *obj)
{
   TObjArray::Add(obj); 
   fCurSitePtr = static_cast<TVKalSite *>(obj);
}
#endif
