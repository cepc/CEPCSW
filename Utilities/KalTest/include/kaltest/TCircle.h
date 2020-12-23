#ifndef TCIRCLE_H
#define TCIRCLE_H
//*************************************************************************
//* ====================
//*  TCircle Class
//* ====================
//*
//* (Description)
//*   A class to implement a circle object.
//* (Requires)
//*     TVCurve
//* (Provides)
//*     class TCircle
//* (Update Recored)
//*   2003/10/03  K.Fujii       Original version.
//*
//*************************************************************************
//
#include <iostream>
#include "TVector2.h"
#include "TVCurve.h"

using namespace std;
//_____________________________________________________________________
//  -----------------------------------
//  Circle Class
//  -----------------------------------

class TCircle : public TVCurve {
public:
   TCircle(Double_t r = 1., Double_t xc = 0., Double_t yc = 0.);
   virtual ~TCircle() {}

   virtual Int_t   CalcXingPointWith(const TCircle  &c,
                                           TVector2  xx[],
                                           Double_t  eps = 1.e-8) const;

   inline virtual       Double_t   GetRadius() const { return fR;  }
   inline virtual const TVector2 & GetCenter() const { return fXc; }

   inline virtual       void       DebugPrint() const;

private:
   Double_t fR;         // radius
   TVector2 fXc;        // center
 
   ClassDef(TCircle,1)      // circle class
};

void TCircle::DebugPrint() const
{
   cerr << " radius = " << fR
        << " center = (" << fXc.X() << "," << fXc.Y() << ")" << endl;
}

#endif
