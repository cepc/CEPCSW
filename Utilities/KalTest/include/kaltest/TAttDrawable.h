#ifndef TATTDRAWABLE_H
#define TATTDRAWABLE_H
//*************************************************************************
//* ===================
//*  TAttDrawable Class
//* ===================
//*
//* (Description)
//*    TAttDrawable class adds drawable attribute to an object.
//* (Requires)
//* 	none
//* (Provides)
//* 	class TAttDrawable
//* (Update Recored)
//*    2004/11/04  K.Fujii	Original very primitive version.
//*
//*************************************************************************
//
#include <Rtypes.h>
//_____________________________________________________________________
//  ------------------------------
//  Base Class for Drawale Objects
//  ------------------------------
//
class TAttDrawable {
public:
   TAttDrawable() {}
   virtual ~TAttDrawable() {}

   virtual void   Draw(const Char_t *opt="");
   virtual void   Draw(Int_t /* color */, const Char_t *opt="");
private:

   ClassDef(TAttDrawable, 1)  // Base class for drawable objects
};

#endif
