#ifndef TVCURVE_H
#define TVCURVE_H
//*************************************************************************
//* ====================
//*  TVCurve Class
//* ====================
//*
//* (Description)
//*   This is the base class for various curves.
//* (Requires)
//*     TObject;
//* (Provides)
//*     class TVCurve
//* (Update Recored)
//*   2003/10/03  K.Fujii       Original version.
//*
//*************************************************************************
//
#include "TObject.h"
//_____________________________________________________________________
//  -----------------------------------
//  Base Class for any curve
//  -----------------------------------

class TVCurve : public TObject {
public:
private:
 
   ClassDef(TVCurve,1)      // Base class for any curve
};

#endif
