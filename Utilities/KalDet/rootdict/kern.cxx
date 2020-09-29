// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME dIcefsdIhiggsdIfucddIKey4hepdICEPCSWdIUtilitiesdIKalDetdIrootdictdIkern
#define R__NO_DEPRECATION

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// The generated code does not explicitly qualifies STL entities
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "/cefs/higgs/fucd/Key4hep/CEPCSW/Utilities/KalDet/src/kern/EXVKalDetector.h"
#include "/cefs/higgs/fucd/Key4hep/CEPCSW/Utilities/KalDet/src/kern/EXVMeasLayer.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void delete_EXVKalDetector(void *p);
   static void deleteArray_EXVKalDetector(void *p);
   static void destruct_EXVKalDetector(void *p);
   static Long64_t merge_EXVKalDetector(void *obj, TCollection *coll,TFileMergeInfo *info);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::EXVKalDetector*)
   {
      ::EXVKalDetector *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::EXVKalDetector >(0);
      static ::ROOT::TGenericClassInfo 
         instance("EXVKalDetector", ::EXVKalDetector::Class_Version(), "", 43,
                  typeid(::EXVKalDetector), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::EXVKalDetector::Dictionary, isa_proxy, 4,
                  sizeof(::EXVKalDetector) );
      instance.SetDelete(&delete_EXVKalDetector);
      instance.SetDeleteArray(&deleteArray_EXVKalDetector);
      instance.SetDestructor(&destruct_EXVKalDetector);
      instance.SetMerge(&merge_EXVKalDetector);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::EXVKalDetector*)
   {
      return GenerateInitInstanceLocal((::EXVKalDetector*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::EXVKalDetector*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void delete_EXVMeasLayer(void *p);
   static void deleteArray_EXVMeasLayer(void *p);
   static void destruct_EXVMeasLayer(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::EXVMeasLayer*)
   {
      ::EXVMeasLayer *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::EXVMeasLayer >(0);
      static ::ROOT::TGenericClassInfo 
         instance("EXVMeasLayer", ::EXVMeasLayer::Class_Version(), "", 106,
                  typeid(::EXVMeasLayer), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::EXVMeasLayer::Dictionary, isa_proxy, 4,
                  sizeof(::EXVMeasLayer) );
      instance.SetDelete(&delete_EXVMeasLayer);
      instance.SetDeleteArray(&deleteArray_EXVMeasLayer);
      instance.SetDestructor(&destruct_EXVMeasLayer);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::EXVMeasLayer*)
   {
      return GenerateInitInstanceLocal((::EXVMeasLayer*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::EXVMeasLayer*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr EXVKalDetector::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *EXVKalDetector::Class_Name()
{
   return "EXVKalDetector";
}

//______________________________________________________________________________
const char *EXVKalDetector::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::EXVKalDetector*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int EXVKalDetector::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::EXVKalDetector*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *EXVKalDetector::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::EXVKalDetector*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *EXVKalDetector::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::EXVKalDetector*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr EXVMeasLayer::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *EXVMeasLayer::Class_Name()
{
   return "EXVMeasLayer";
}

//______________________________________________________________________________
const char *EXVMeasLayer::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::EXVMeasLayer*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int EXVMeasLayer::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::EXVMeasLayer*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *EXVMeasLayer::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::EXVMeasLayer*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *EXVMeasLayer::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::EXVMeasLayer*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void EXVKalDetector::Streamer(TBuffer &R__b)
{
   // Stream an object of class EXVKalDetector.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(EXVKalDetector::Class(),this);
   } else {
      R__b.WriteClassBuffer(EXVKalDetector::Class(),this);
   }
}

namespace ROOT {
   // Wrapper around operator delete
   static void delete_EXVKalDetector(void *p) {
      delete ((::EXVKalDetector*)p);
   }
   static void deleteArray_EXVKalDetector(void *p) {
      delete [] ((::EXVKalDetector*)p);
   }
   static void destruct_EXVKalDetector(void *p) {
      typedef ::EXVKalDetector current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around the merge function.
   static Long64_t  merge_EXVKalDetector(void *obj,TCollection *coll,TFileMergeInfo *) {
      return ((::EXVKalDetector*)obj)->Merge(coll);
   }
} // end of namespace ROOT for class ::EXVKalDetector

//______________________________________________________________________________
void EXVMeasLayer::Streamer(TBuffer &R__b)
{
   // Stream an object of class EXVMeasLayer.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(EXVMeasLayer::Class(),this);
   } else {
      R__b.WriteClassBuffer(EXVMeasLayer::Class(),this);
   }
}

namespace ROOT {
   // Wrapper around operator delete
   static void delete_EXVMeasLayer(void *p) {
      delete ((::EXVMeasLayer*)p);
   }
   static void deleteArray_EXVMeasLayer(void *p) {
      delete [] ((::EXVMeasLayer*)p);
   }
   static void destruct_EXVMeasLayer(void *p) {
      typedef ::EXVMeasLayer current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::EXVMeasLayer

namespace {
  void TriggerDictionaryInitialization_kern_Impl() {
    static const char* headers[] = {
"0",
0
    };
    static const char* includePaths[] = {
"/cefs/higgs/fucd/Key4hep/CEPCSW/Utilities/KalTest",
"/cefs/higgs/fucd/Key4hep/CEPCSW/Utilities/KalDet/src/gen",
"/cefs/higgs/fucd/Key4hep/CEPCSW/Utilities/KalDet/src/kern",
"/cefs/higgs/fucd/Key4hep/CEPCSW/Utilities/KalDet/src/lctpc/gearTPC",
"/cvmfs/sft.cern.ch/lcg/releases/ROOT/v6.20.02-d9e99/x86_64-slc6-gcc8-opt/include/",
"/cefs/higgs/fucd/Key4hep/CEPCSW/build/Utilities/KalDet/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "kern dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
class __attribute__((annotate(R"ATTRDUMP(Abstract measurement layer class)ATTRDUMP"))) EXVKalDetector;
class __attribute__((annotate(R"ATTRDUMP(Abstract measurement layer class)ATTRDUMP"))) EXVMeasLayer;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "kern dictionary payload"

#ifndef HANDLE_DICT_EXCEPTIONS
  #define HANDLE_DICT_EXCEPTIONS IGNORED_FOR_CINT
#endif

#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
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

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"EXVKalDetector", payloadCode, "@",
"EXVMeasLayer", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("kern",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_kern_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_kern_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_kern() {
  TriggerDictionaryInitialization_kern_Impl();
}
