#ifndef TATTELEMENT_H
#define TATTELEMENT_H
//*************************************************************************
//* ===================
//*  TAttElement Class
//* ===================
//*
//* (Description)
//*    TAttElement class adds constituent attribute to an object.
//* (Requires)
//* 	none
//* (Provides)
//* 	class TAttElement
//* (Update Recored)
//*    2003/10/10  K.Fujii	Original very primitive version.
//*
//*************************************************************************
//
#include <Rtypes.h>

//_____________________________________________________________________
//  --------------------------------
//  Base Class for Element Objects
//  --------------------------------
//
class TAttElement {
public:
   TAttElement() : fParentPtr(0) {}
   virtual ~TAttElement() {}

   inline virtual const TAttElement  & GetParent(Bool_t recur = kTRUE) const;

   inline virtual void  SetParentPtr(TAttElement *obj)  { fParentPtr = obj; }

private:
   TAttElement *fParentPtr;	 // pointer to parent

   ClassDef(TAttElement,1)  // Base class for lockable objects
};

//_____________________________________________________________________
//  --------------------------------
//  Inline functions, if any
//  --------------------------------
const TAttElement & TAttElement::GetParent(Bool_t recursive) const
{
   if (fParentPtr) {
      if (recursive) return fParentPtr->GetParent(recursive);
      else           return *fParentPtr;
   } else {
      return *this;
   }
}

#endif
