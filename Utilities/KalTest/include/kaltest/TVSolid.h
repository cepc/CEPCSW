#ifndef TVSOLID_H
#define TVSOLID_H
//*************************************************************************
//* ====================
//*  TVSolid Class
//* ====================
//*
//* (Description)
//*   This is the base class for various solids.
//* (Requires)
//*     TObject;
//* (Provides)
//*     class TVSolid
//* (Update Recored)
//*   2003/10/03  K.Fujii       Original version.
//*
//*************************************************************************
//
#include "TObject.h"
//_____________________________________________________________________
//  -----------------------------------
//  Base Class for any solid
//  -----------------------------------

class TVSolid : public TObject {
public:
private:
 
   ClassDef(TVSolid,1)      // Base class for any solid
};

#endif
