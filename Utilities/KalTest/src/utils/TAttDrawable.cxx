//*************************************************************************
//* ====================
//*  TAttDrawable Class
//* ====================
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
#include "TAttDrawable.h"
//_________________________________________________________________________
//  -------------------------------
//  Base Class for Drawable Objects
//  -------------------------------
//
ClassImp(TAttDrawable)

//=========================================================================
//* Draw ------------------------------------------------------------------

void TAttDrawable::Draw(const Char_t *opt)
{
   static Int_t color = 0;
   color++;
   color %= 10;
   Draw(color, opt);
}

void TAttDrawable::Draw(Int_t /* color */, const Char_t * /* opt */)
{
}
