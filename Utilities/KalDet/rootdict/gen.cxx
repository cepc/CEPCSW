// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME dIcefsdIhiggsdIfucddIKey4hepdICEPCSWdIUtilitiesdIKalDetdIrootdictdIgen
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
#include "/cefs/higgs/fucd/Key4hep/CEPCSW/Utilities/KalDet/src/gen/EXEventGen.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static TClass *EXEventGen_Dictionary();
   static void EXEventGen_TClassManip(TClass*);
   static void delete_EXEventGen(void *p);
   static void deleteArray_EXEventGen(void *p);
   static void destruct_EXEventGen(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::EXEventGen*)
   {
      ::EXEventGen *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::EXEventGen));
      static ::ROOT::TGenericClassInfo 
         instance("EXEventGen", "", 13,
                  typeid(::EXEventGen), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &EXEventGen_Dictionary, isa_proxy, 4,
                  sizeof(::EXEventGen) );
      instance.SetDelete(&delete_EXEventGen);
      instance.SetDeleteArray(&deleteArray_EXEventGen);
      instance.SetDestructor(&destruct_EXEventGen);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::EXEventGen*)
   {
      return GenerateInitInstanceLocal((::EXEventGen*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::EXEventGen*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *EXEventGen_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::EXEventGen*)0x0)->GetClass();
      EXEventGen_TClassManip(theClass);
   return theClass;
   }

   static void EXEventGen_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrapper around operator delete
   static void delete_EXEventGen(void *p) {
      delete ((::EXEventGen*)p);
   }
   static void deleteArray_EXEventGen(void *p) {
      delete [] ((::EXEventGen*)p);
   }
   static void destruct_EXEventGen(void *p) {
      typedef ::EXEventGen current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::EXEventGen

namespace {
  void TriggerDictionaryInitialization_gen_Impl() {
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
#line 1 "gen dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
class EXEventGen;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "gen dictionary payload"

#ifndef HANDLE_DICT_EXCEPTIONS
  #define HANDLE_DICT_EXCEPTIONS IGNORED_FOR_CINT
#endif

#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#ifndef __EXEVENTGEN__
#define __EXEVENTGEN__

#include "kaltest/TKalDetCradle.h"
#include "kaltest/THelicalTrack.h"
#include "TMath.h"

class EXEventGen {
public:
   EXEventGen(TKalDetCradle const &cradle, TObjArray &kalhits)
             : fCradlePtr(&cradle), fHitBufPtr(&kalhits) {}
   virtual ~EXEventGen() {}

   THelicalTrack GenerateHelix(Double_t pt,
                               Double_t cosmin,
                               Double_t cosmax,
                               Double_t phimin=0.,
                               Double_t phimax=2*TMath::Pi(),
                               TVector3 xv0=TVector3(0.,0.,0.));
   void          Swim(THelicalTrack &heltrk);

   static void     SetT0(Double_t t0) { fgT0 = t0;   }
   static Double_t GetT0()            { return fgT0; }

private:
   TKalDetCradle const *fCradlePtr;     // pointer to detector system
   TObjArray     *fHitBufPtr;     // pointer to hit array

   static Double_t  fgT0;         // t0

};

#endif

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"EXEventGen", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("gen",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_gen_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_gen_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_gen() {
  TriggerDictionaryInitialization_gen_Impl();
}
