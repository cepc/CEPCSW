#ifndef TKALFILTERCOND_H
#define TKALFILTERCOND_H
//*************************************************************************
//* =====================
//*  TKalFilterCond Class
//* =====================
//*
//* (Description)
//*   A class to specify filter conditions used in Kalman filter.
//* (Requires)
//* (Provides)
//* 	class TKalFilterCond
//* (Update Recored)
//*   2010/04/06  K.Fujii        Original Version.
//*
//*************************************************************************
#include "Rtypes.h"
//_____________________________________________________________________
//  ------------------------------
//   Filter condition class
//  ------------------------------

class TKalTrackSite;

class TKalFilterCond {
 public:
  
  // need virtual destructor is we have virtual functions
  virtual ~TKalFilterCond() {};
  
  virtual Bool_t IsAccepted(const TKalTrackSite &site);
  
  ClassDef(TKalFilterCond,1)  // Base class for detector system
    };
#endif
