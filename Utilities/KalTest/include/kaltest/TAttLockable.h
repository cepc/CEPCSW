#ifndef TATTLOCKABLE_H
#define TATTLOCKABLE_H
//*************************************************************************
//* ===================
//*  TAttLockable Class
//* ===================
//*
//* (Description)
//*    TAttLockable class adds lockable attribute to an object.
//* (Requires)
//* 	none
//* (Provides)
//* 	class Lockable
//* (Update Recored)
//*    1999/06/05  K.Fujii	Original very primitive version.
//*
//*************************************************************************
//
#include <Rtypes.h>
//_____________________________________________________________________
//  ------------------------------
//  Base Class for Lockale Objects
//  ------------------------------
//
class TAttLockable {
public:
   TAttLockable() : fStatus(kFALSE) {}
   virtual ~TAttLockable() {}

   inline virtual Bool_t IsLocked() const { return fStatus;   }
   inline virtual void   Lock()           { fStatus = kTRUE;  }
   inline virtual void   Unlock()         { fStatus = kFALSE; }
private:
   Bool_t fStatus;	 // lock byte

   ClassDef(TAttLockable,1)  // Base class for lockable objects
};

#endif
