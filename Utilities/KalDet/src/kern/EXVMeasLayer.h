#ifndef EXVMEASLAYER_H
#define EXVMEASLAYER_H
//*************************************************************************
//* ====================
//*  EXVMeasLayer Class
//* ====================
//*
//* (Description)
//*   Abstract measurement layer class used by TVTrackHit
//* (Requires)
//*     TVMeasLayer
//* (Provides)
//*     class EXVMeasLayer
//* (Update Recored)
//*   2009/11/23  K.Ikematsu   Derived from KalTest/examples/kaltest/
//*                                         hybrid/kern/EXVMeasLayer.h
//*
//* $Id: EXVMeasLayer.h,v 1.1.1.1 2009-11-24 00:13:59 ikematsu Exp $
//*************************************************************************
//
#include <TVector3.h>
#include <kaltest/TKalMatrix.h>
#include <kaltest/TCylinder.h>
#include <kaltest/TVMeasLayer.h>
#include <kaltest/TAttDrawable.h>
#include <kaltest/KalTrackDim.h>
#include <TString.h>

class TVTrackHit;
#include <TNode.h>

class EXVMeasLayer : public TVMeasLayer, public TAttDrawable {

public:
  static Bool_t kActive;
  static Bool_t kDummy;

  // Ctors and Dtor

  EXVMeasLayer(TMaterial &min,
               TMaterial &mout,
               Bool_t type = EXVMeasLayer::kActive,
               const Char_t *name = "MeasL");
  virtual ~EXVMeasLayer();

  virtual void ProcessHit(const TVector3 &xx,
                          TObjArray &hits) const = 0;

  inline TString  GetMLName () const { return fName;    }
  inline TNode   *GetNodePtr() const { return fNodePtr; }

  inline void     SetNodePtr(TNode *nodep) { fNodePtr = nodep; }

private:
  TString  fName;     // layer name
  TNode   *fNodePtr;  // node pointer

  ClassDef(EXVMeasLayer, 1)  // Abstract measurement layer class
};
#endif
