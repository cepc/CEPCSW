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
//*                                         hybrid/kern/EXVMeasLayer.cxx
//*
//* $Id: EXVMeasLayer.cxx,v 1.1.1.1 2009-11-24 00:13:59 ikematsu Exp $
//*************************************************************************
//
#include "EXVMeasLayer.h"
#include <TNode.h>

Bool_t EXVMeasLayer::kActive = kTRUE;
Bool_t EXVMeasLayer::kDummy  = kFALSE;

ClassImp(EXVMeasLayer)

EXVMeasLayer::EXVMeasLayer(TMaterial &min,
                           TMaterial &mout,
                           Bool_t     isactive,
                     const Char_t    *name)
            : TVMeasLayer(min, mout, isactive),
              fName(name),
              fNodePtr(0)
{
}

EXVMeasLayer::~EXVMeasLayer()
{
}
