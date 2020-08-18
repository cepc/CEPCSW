#ifndef EXVKALDETECTOR_H
#define EXVKALDETECTOR_H
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
//*                                         hybrid/kern/EXVKalDetector.h
//*
//* $Id: EXVKalDetector.h,v 1.1.1.1 2009-11-24 00:13:59 ikematsu Exp $
//*************************************************************************
//
#include <TVector3.h>
#include <kaltest/TVKalDetector.h>
#include <kaltest/TAttDrawable.h>

class TNode;

/**
 * Base class to make a detector drawable, add a magnetic field,
 * a power switch (whatever the use may be).
 * 
 * Killenb: I removed the TAttDrawable for the moment. The TNode pointer 
 * stuff and the implementation of Draw belong to the TAttDrawable anyway. So if 
 * the drawability is needed move it to TAttDrawable and just inherit from it.
 *
 * \deprecated EXVKalDetector
 */
//class EXVKalDetector : public TVKalDetector, public TAttDrawable {
class EXVKalDetector : public TVKalDetector {

public:
  EXVKalDetector(Double_t bField, Int_t m = 100);
  virtual ~EXVKalDetector();

  /// Return whether the power is on. Currently hard coded to true.
  inline virtual Bool_t IsPowerOn() const { return true;   }

  /// Turn the power on. Currently ignored.
  inline virtual void   PowerOn  ()       { fIsPowerOn = kTRUE;  }

  /// Turn the power off. Currently ignored.
  inline virtual void   PowerOff ()       { fIsPowerOn = kFALSE; }

  /// Returns a single double value with a 3D point as an input. 
  /// Completely unphysical interface. Either the magnetic field varies with the position,
  /// in which case you need a three-dimensional return value, or B can be desrcibed as single 
  /// value, which means it is homogeneous and thus indepenent from the position.
  /// Currently it does the only reasonable thing: It ignores the argument and returns the 
  /// constant value given in the constructor.
  virtual Double_t GetBfield (const TVector3 &xx = TVector3(0.,0.,0.)) const
                                          { return fBfield; }

protected:
  Bool_t fIsPowerOn;           // power status
  Double_t  fBfield;   // magnetic field [T]

  ClassDef(EXVKalDetector, 1)  // Abstract measurement layer class
};
#endif
