//*************************************************************************
//* ======================
//*  EXVKalDetector Class
//* ======================
//*
//* (Description)
//*   Abstract detector class for Kalman filter
//* (Requires)
//*     TVKalDetector
//* (Provides)
//*     class EXVKalDetector
//* (Update Recored)
//*   2009/11/23  K.Ikematsu   Derived from KalTest/examples/kaltest/
//*                                         hybrid/kern/EXVKalDetector.cxx
//*   2010/11/17  K.Fujii      Changed unit system to (mm, nsec, T) 
//*                            from (cm, nsec, kG)
//*
//* $Id: EXVKalDetector.cxx,v 1.1.1.1 2009-11-24 00:13:59 ikematsu Exp $
//*************************************************************************
//
#include "EXVKalDetector.h"

ClassImp(EXVKalDetector)

EXVKalDetector::EXVKalDetector(Double_t bField, Int_t m)
  : TVKalDetector(m),
    fIsPowerOn(kTRUE),
    fBfield(bField)
{
}

EXVKalDetector::~EXVKalDetector()
{
}
