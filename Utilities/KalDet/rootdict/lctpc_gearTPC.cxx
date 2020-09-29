// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME dIcefsdIhiggsdIfucddIKey4hepdICEPCSWdIUtilitiesdIKalDetdIrootdictdIlctpc_gearTPC
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
#include "/cefs/higgs/fucd/Key4hep/CEPCSW/Utilities/KalDet/src/lctpc/gearTPC/EXTPCHit.h"
#include "/cefs/higgs/fucd/Key4hep/CEPCSW/Utilities/KalDet/src/lctpc/gearTPC/EXTPCKalDetector.h"
#include "/cefs/higgs/fucd/Key4hep/CEPCSW/Utilities/KalDet/src/lctpc/gearTPC/EXTPCMeasLayer.h"
#include "/cefs/higgs/fucd/Key4hep/CEPCSW/Utilities/KalDet/src/lctpc/gearTPC/GearTPCCylinderHit.h"
#include "/cefs/higgs/fucd/Key4hep/CEPCSW/Utilities/KalDet/src/lctpc/gearTPC/GearTPCCylinderMeasLayer.h"
#include "/cefs/higgs/fucd/Key4hep/CEPCSW/Utilities/KalDet/src/lctpc/gearTPC/GearTPCHit.h"
#include "/cefs/higgs/fucd/Key4hep/CEPCSW/Utilities/KalDet/src/lctpc/gearTPC/GearTPCKalDetector.h"
#include "/cefs/higgs/fucd/Key4hep/CEPCSW/Utilities/KalDet/src/lctpc/gearTPC/GearTPCMeasLayer.h"

// Header files passed via #pragma extra_include

namespace {
  void TriggerDictionaryInitialization_lctpc_gearTPC_Impl() {
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
#line 1 "lctpc_gearTPC dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "lctpc_gearTPC dictionary payload"

#ifndef HANDLE_DICT_EXCEPTIONS
  #define HANDLE_DICT_EXCEPTIONS IGNORED_FOR_CINT
#endif

#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#ifndef LCTPC_EXTPCHIT_H
#define LCTPC_EXTPCHIT_H

#include "GearTPCCylinderHit.h"
#include <kaltest/TVMeasLayer.h>

/**
 * A backward compatibility class for GearTPCCylinderHit.
 * Do not use this in new code, but use GearTPCCylinderHit directly. 
 * This class extends the GearTPCCylinderHit by a side, which is never used anywhere.
 *
 * \deprecated EXTPCHit
 */

class EXTPCHit : public kaldet::GearTPCCylinderHit
{
  public:
  /// The default constructor. 
  EXTPCHit(Int_t m = kMdim);

  /// Constructor initialising the original hit as 3D coordinates.
  EXTPCHit(const TVMeasLayer &ms,
                 Double_t       *x,
                 Double_t       *dx,
                 Int_t           side,
                 Double_t        v,
           const TVector3       &xx,
                 Double_t        b,
	         Int_t           m = kMdim);

  /// Constructor initialising the original hit with a reference pointer.
  EXTPCHit(const TVMeasLayer &ms,
                 Double_t       *x,
                 Double_t       *dx,
                 Int_t           side,
                 Double_t        v,
           const void           *hitp,
                 Double_t        b,
                 Int_t           m = kMdim);

  /// The destructor.
  virtual ~EXTPCHit();

  /// Get the side value which has been set in the constructor.
  inline       Int_t     GetSide  () const { return fSide;   }

 private: 
  Int_t           fSide;    /// (-1, +1) = (-z side, +z side)

  //  ClassDef(EXTPCHit, 1)  // EXTPC hit class

};

#endif // LCTPC_EXTPCHIT_H
#ifndef LCTPC_EXTPCKALDETECTOR_H
#define LCTPC_EXTPCKALDETECTOR_H

#include "GearTPCKalDetector.h"

/**
 *  A backward compatibility class for GearTPCKalDetector.
 *  It basically provides a static instance of the detector which can be
 *  accessed via the GetInstance() method.
 *  In addition it provides the static GetVDrift() and GetBField(), which are used
 *  in the old code. The use of this class is highly depreciated.
 *
 *  \deprecated EXTPCKalDetector
 */

class EXTPCKalDetector: public kaldet::GearTPCKalDetector
{
  private:
  /// As this this a singleton class the constructor is private.
  EXTPCKalDetector(const gear::GearMgr& gearMgr);

public:
  /// The destructor.
  virtual ~EXTPCKalDetector();
  
  /// Static access function to the singleton instance.
  static EXTPCKalDetector * GetInstance();

  /// Returns the hard coded drift velocity of 76.e-3 mm/ns. 
  static Double_t GetVdrift() { return fgVdrift; }

  /// Static function to access the magnetic field.
  static Double_t GetBfield();

 private:
  static Double_t           fgVdrift;   //< The drift velocity.
  static EXTPCKalDetector * fgInstance; //< The singleton pointer.

  Double_t fBField; //< The magnetic field

  //  ClassDef(EXTPCKalDetector, 1)  // User defined detector class

};

#endif // LCTPC_EXTPCKALDETECTOR_H
#ifndef LCTPC_EXTPCMEASLAYER_H
#define LCTPC_EXTPCMEASLAYER_H

#include "GearTPCMeasLayer.h"
#include <kaltest/TCylinder.h>

/**
 *  A backward compatibility class for GearTPCCylinderMeasLayer.
 *  It introduces one module and one row, which are associated to the layer. 
 *  This is deprecated, the GearTPCCylinderMeasLayer provides an assiciation of several
 *  module-row pairs to the layer.
 
 *  This class is an intermediate inheritance class so a GearTPCCylinderMeasLayer can be 
 *  instantiated (as should be in the current code) and the old code can cast the 
 *  TObject pointer, which is delivered by the detector cradle, to an EXTPCMeasLayer.
 *
 *  \attention Do not use any of these function in new code. All new code should still run 
 *  after this class has been removed from the ineritance chain.
 *
 * \deprecated EXTPCMeasLayer
 */
class EXTPCMeasLayer : public kaldet::GearTPCMeasLayer, public TCylinder
{

public:
  /// Minimal constructor for this (partially) virtual class.
  EXTPCMeasLayer(TMaterial &min,
		 TMaterial &mout,
		 Int_t      module,
		 Int_t      row,
		 Double_t   r0,
		 Double_t   lhalf,
		 TVector3   xc,
		 Bool_t     isPerfect,
		 Bool_t     isActive,
		 Double_t   sigmaX0 = 0., //< the constant part of sigmaX
		 Double_t   sigmaX1 = 0., //< the z-dependent part of sigmaX
		 Double_t   sigmaZ0 = 0., //< the constant part of sigmaZ
		 Double_t   sigmaZ1 = 0.); //< the z-dependent part of sigmaZ

  /// The destructor.
  virtual ~EXTPCMeasLayer();

  /**
   *  Get the module associated with this layer (deprecated).
   *  \attention Do not programme against this when using the GearTPC interface. 
   *  This is for  backward compatibility only!!!
   */
  Int_t GetModuleID() const;
  
  /** 
   *  Get the layer ID (i.\ e.\ row in the module) associated with this Kalman layer (deprecated).
   *
   *  \attention Do not programme against this when using the GearTPC interface. 
   * This is for  backward compatibility only!!!
   */
  Int_t GetLayerID () const;

  /** Deprecated XvToMv which in addition to the position takes a side. 
   *  Side is ignored and XvToMv without the side is called.
   * \attention Do not programme against this when using the GearTPC interface. 
   * This is for  backward compatibility only!!!
   */
  TKalMatrix XvToMv(const TVector3 &xv, Int_t side) const;

  /** The fully virtual declaration of XvToMv. It is called within the version which also takes
   *  the side argument, but is implemented in GearTPCCylinderMeasLayer.
   */
  virtual TKalMatrix XvToMv(const TVector3 &xv) const = 0;
  
  /** Smear the incoming hit in the layes measurement surface and place the result into the TObjArray
   *  which is given as argument.
   *  From a design point of view this function should not be in the detector class but in a 
   *  simulation  extension. It is only put in for compatibility reasons.
   *  \attention Do not programme against this when using the GearTPC interface. 
   * This is for  backward compatibility only!!!
   */
  virtual void ProcessHit(const TVector3  &xx, TObjArray &hits) const;

  //ClassDef(EXTPCMeasLayer, 1)  // User defined measurement layer class
};

#endif // LCTPC_EXTPCMEASLAYER_H
#ifndef GEARTPCCYLINDERHIT_H
#define GEARTPCCYLINDERHIT_H

#include <kaltest/KalTrackDim.h>
#include "GearTPCHit.h"
#include <kaltest/TVMeasLayer.h>

namespace kaldet{

/** The cylindrical implementation of the GearTPCHit.
 */
class GearTPCCylinderHit : public GearTPCHit {

public:
    /// KILLENB What does this constructor do? Best throw it out, it does not 
    /// properly initialise the class at all, does it?
  GearTPCCylinderHit(Int_t m = kMdim);

  /** Constructor to initialise the GearTPCHit using space point coordinates (TVector3) as original hit.
   */
  GearTPCCylinderHit(const TVMeasLayer &ms,
                 Double_t       *x,
                 Double_t       *dx,
           const TVector3       &xx,
                 Double_t        b,
                 Double_t        v,
                 Int_t           m = kMdim);

  /** Constructor using a pointer to the original hit as reference.
   */
  GearTPCCylinderHit(const TVMeasLayer &ms,
                 Double_t       *x,
                 Double_t       *dx,
           const void           *hitp,
                 Double_t        b,
                 Double_t        v,
                 Int_t           m = kMdim);

  /** The dectructor.
   */
  virtual ~GearTPCCylinderHit();

  /** Implementation of the space vector (xv) to measurement vector (mv) calculation
   *  for a cylindrical hit.
   */
  virtual TKalMatrix XvToMv(const TVector3 &xv, Double_t t0) const;
  
  /** Print some debug output to std err.
   */
  virtual void       DebugPrint(Option_t *opt = "")          const;
};

}//namespace kaldet

#endif //GEARTPCCYLINDERHIT_H
#ifndef GEARTPCCYLINDERMEASLAYER_H
#define GEARTPCCYLINDERMEASLAYER_H
#include <TVector3.h>
#include <kaltest/TKalMatrix.h>
#include <kaltest/TCylinder.h>
#include <EXTPCMeasLayer.h>
//#include <KalTrackDim.h>

#include <TMath.h>

#include <set>

class TVTrackHit;

namespace kaldet
{

  /**
   *  A cylindrical measurement layer.
   */
  class GearTPCCylinderMeasLayer 
    : public EXTPCMeasLayer
    /* this is the original code which should be reactivated once the EXTPCMeasLayer is phased out:
    : public GearTPCMeasLayer, public TCylinder 
    */
  {
    
  public:
    /** The constructor.
     *  If the layer is perfect it is always a full circle. The constructor forces 
     *  phiMin and phiMax to +-TMath::Pi(). Differing values will be ignored and a warning is
     *  printed.
     *
     *  Note: The current implementation forces the layer to be perfect. Segmented layers are
     *  not supported yet. 
     *
     *  Note: for backward compatibility this is derrived from EXTPCMeasLayer.
     *  After the change to the GearTPC interface this should be changed to GearTPCMeasLayer and
     *  TCylinder, as EXTPCMeasLayer is inherrited now.
     *  The current version ensures compatibility for the transition phase.
     */ 
  GearTPCCylinderMeasLayer(TMaterial &min,
			   TMaterial &mout,
			   Int_t      module,
			   Int_t      row,
			   Double_t   r0,
			   Double_t   lhalf,
			   TVector3   xc =  TVector3(),
			   Bool_t     isPerfect = true,
			   Bool_t     isActive = true,
			   Double_t   sigmaX0 = 0.,
			   Double_t   sigmaX1 = 0.,
			   Double_t   sigmaZ0 = 0.,
			   Double_t   sigmaZ1 = 0.,
			   Double_t   phiMin = -TMath::Pi(),
			   Double_t   phiMax = TMath::Pi());

  /**
   * The desctructor.
   */
  virtual ~GearTPCCylinderMeasLayer();

  // Parent's pure virtuals that must be implemented

  /** Implements kaltest::TVMeasLayer's XvToMv. I have no idea why there are two arguments.
   *  It ignores ht and just calls  XvToMv(xv).
   */
  virtual TKalMatrix XvToMv    (const TVTrackHit &ht,
                                const TVector3   &xv)   const;

  /** Implements the coordinate transformation from the space vector xv to the
   *  measurement vector (Kalman matrix).
   */
  virtual TKalMatrix XvToMv    (const TVector3   &xv)   const;

  /** Implements the conversion from a Kalman hit (measurement vector) to 
   *  a 3D space point.
   */
  virtual TVector3   HitToXv   (const TVTrackHit &ht)   const;

  /**
   * Implements CalcDhDa, whatever that is.
   */
  virtual void       CalcDhDa  (const TVTrackHit &ht,
                                const TVector3   &xv,
                                const TKalMatrix &dxphiada,
                                      TKalMatrix &H)    const;
  /** Implements the sorting policy.
   *  The layers are first sorted by radius + offset. This offset is only
   *  useful for segments of a cylinder, like the LP1.
   *  As offsets in this case can be positive or negative, but only make sense in one 
   *  direction (you need a continuous number), we only allow offsets in x.
   *  This should not be too much of a problem, you should be able to rotate your coordinates
   *  so the offset is in x. If not you have to extend the sorting policy. (Please thake
   *  care not to reduce versatility when doing so. You might want to implement your own class?)
   *  
   *  For equal radii  + offset the layers are sorted by moduleID. As we have to squeeze this 
   *  information into only one number, we multiply the radius + offset by 1e9 and add the moduleID.
   *  A double has a precision of 53 bits, which is 15.9 digits. So the radius can be up to 1e6.9 mm
   *  without causing the last digit of the the ModuleID to be cut, and for up to 1000 modules the
   *  layers can be distinguished down to 1 nm without the two numbers mixing, or down to 1 micron
   *  with up to 1.000.000 modules.
   * 
   *  The additional sorting by module is intended for cylinder segments. Here only one module/row
   *  per layer is allowed, so we just take the first entry in the set. In case of a perfect layer
   *  it does not matter because there should only be one layer at this radius, so the sort order
   *  should not be affected by adding an arbitrary module ID (as long as the module ID is < 1e6, as 
   *  described above).
   */
  virtual Double_t   GetSortingPolicy() const;

  /**
    * Creates a GearTPCCylinderHit and hands over the ownership. 
   */
  virtual GearTPCHit * createHit(Double_t * meas,
				 Double_t * dmeas,
				 void * hitPointer, 
				 Double_t bField,
				 Double_t vDrift,
				 Int_t           m = kMdim) const;



protected:
  Double_t fPhiMin;   //< Minimum phi.
  Double_t fPhiMax;   //< Maximum phi.

};

}//namespace kaldet
#endif
#ifndef GEARTPCHIT_H
#define GEARTPCHIT_H

#include <kaltest/KalTrackDim.h>
#include <kaltest/TVTrackHit.h>
#include <kaltest/TVMeasLayer.h>

namespace kaldet{

/** Base class of a hit for GearTPCKalDetector. It extends the TVTrackHit with the functionality to
 *  store a space point or a pointer to the original hit for reference. In addition it stores
 *  the local drift velocity and allows sorting of the hits (according to the distance to the 
 *  origin).
 *
 *  It does not implement the purely virtual functions of the TVTrackHit, which happens in the 
 *  specialisations for cylindrical and planar measurement layers.
 */
class GearTPCHit : public TVTrackHit {

public:
    /// KILLENB What does this constructor do? Best throw it out, it does not 
    /// properly initialise the class at all, does it?
  GearTPCHit(Int_t m = kMdim);

  /** Constructor to initialise the GearTPCHit using space point coordinates (TVector3) as original hit.
   */
  GearTPCHit(const TVMeasLayer &ms,
                 Double_t       *x,
                 Double_t       *dx,
           const TVector3       &xx,
                 Double_t        b,
                 Double_t        v,
                 Int_t           m = kMdim);

  /** Constructor using a pointer to the original hit as reference.
   */
  GearTPCHit(const TVMeasLayer &ms,
                 Double_t       *x,
                 Double_t       *dx,
           const void           *hitp,
                 Double_t        b,
                 Double_t        v,
                 Int_t           m = kMdim);

  /** The dectructor.
   */
  virtual ~GearTPCHit();

  /**
   * The sorting policy of hits is implemented as the distance to the origin.
   *
   * Note: The sorting of hits does not necessarily correspond to the sort order of 
   * the corresponding Kalman layer!
   */
  virtual Double_t   GetSortingPolicy()                      const;
  
  /**
   * Compare two hits according to their sorting policy.
   * Returns
   * \li -1 if this hits sorting policy is smaller
   * \li 0  if both soting policies are equal
   * \li 1 if this hits hits sortin policy is larger
   *
   * Killenb: n.b. Who comes up with this wierd stuff? Ever head of anything like a 
   * `less than operator` or `comparison operator`?
   */
  virtual Int_t      Compare(const TObject *obj)             const;

  /**
   * Returns true.
   */
  virtual Bool_t     IsSortable()                            const;
  
  /// Get the pointer to the reference hit. 0 if the TVector3 has been used for initialisation.
  inline const void     *GetHitPtr() const { return fHitPtr; }

  /// Get the referece position. (0, 0, 0)  if the reference pointer has been used for initialisation.
  inline       TVector3  GetExactX() const { return *fXXPtr; }

  /// Get the local drift velocity set in the constructor.
  inline       Double_t  GetVdrift() const { return fVDrift; }

protected:
  const TVector3 *fXXPtr;   //< pointer to exact hit
  const void     *fHitPtr;  //< pointer to raw Hit object

  Double_t        fVDrift;  //< the local drift veclocity at this point
  

};

}//namespace kaldet

#endif
#ifndef GEARTPCKALDETECTOR_H
#define GEARTPCKALDETECTOR_H

#include "kaltest/TVKalDetector.h"

#include "GearTPCMeasLayer.h"

#include <map>

namespace gear{
  class GearMgr ;
}

namespace kaldet{

  /**
   * The LCTPC implementation for a TPC which is completely instantiated from GEAR.
   * 
   */
class GearTPCKalDetector : public TVKalDetector {

public:
    /** 
     * The constructor. All information to initialise the TPC is taken from GEAR.
     *
     * As a pragmatic approach to avoid dealing with conditions data and material databases,
     * the information about the material budget and the resolution of the layers
     * is taken from the GEAR file as user parameters. If the parameters are not found in the
     * file the previously hard coded parameters are used as default, which ensures backward
     * compatibility.
     *
     * The gas properties for the matrial budget can be given as user parameters 
     * for the TPCParameters:
     * \param  TPCGas_A The mean atomic mass (default 36.2740552)
     * \param  TPCGas_Z The mean number of protons (default 16.4)
     * \param  TPCGas_density The density (default 0.749e-3 in which units?)
     * \param  TPCGas_radlen The radiation length (default 2.392e4 in which units?)
     *
     * The default gas parameters (are supposed to) correspond to Ar/CH4 90/10.
     * N.B.: KILLENB: I think there is a bug in the calculation, the mean A should be
     * 37.6 instead of 36.3 (see source code).
     * In addition the description as a single TMaterial is not good. 
     * Using TMixture would be better.
     *
     * The reslution is calculated as \f$\sigma_x = \sqrt{x_0^2 + x_1^2 \cdot z}\f$.
     * This requires z to be proportional to the drift distance, i.\ e. z=0 is at the readout.

     * The resolution of the layers can be given as user parameters in each TPCModule 
     * section of the GEAR xml file.
     * \param sigmax0 The constant part of the x resolution (default 38.3e-3 mm)
     * \param sigmax1 The drift distance dependent part of the x resolution 
     *                (default 6.74e-3 mm/sqrt(mm) )
     * \param sigmaz0 The constant part of the z resolution (default 0.5 mm)
     * \param sigmaz1 The drift distance dependent part the z resolution
     *                (default 10.2e-3 mm/sqrt(mm) )
     */
    GearTPCKalDetector(const gear::GearMgr& gearMgr);

    /// The destructor.
    virtual ~GearTPCKalDetector();

    /**
     * Get access to the measurement layers using moduleID and row.
     * Do not directly access the measurement layers using At() 
     * because the order depends on the order in the gear file.
     * Throws a gear::Exception if the row on the module is not defined.
     */
    virtual GearTPCMeasLayer const * GetMeasLayer(int moduleID, int row) const;

protected:
    /// Map which contains the information which measurement layer is stored
    /// at which position in the array.
    std::map< std::pair<int, int >, Int_t > moduleRowToMeasurementLayerMap;
};

}// namespace kaldet
#endif //GEARTPCKALDETECTOR_H
#ifndef GEARTPC_MEASLAYER_H
#define GEARTPC_MEASLAYER_H

#include <kaltest/TVMeasLayer.h>
#include <set>

namespace kaldet
{

  class GearTPCHit;

  /**
   * The GearTPCMeasLayer class introduces the z-dependent resolutions sigmaX and sigmaZ
   * as well as Gear modules and rows which correspond to this layer.
   *
   * If the layer is defined as a perfect layer this means all modules are perfectly alligned 
   * and more than one module/row can be assigned to this layer. You can add them using AddModuleRow.
   * The perfect layer should contain all the moduleRows on it, so it is guaranteed that the
   * user can access all neighbouring modules this way. 
   *
   * If the layer is not defined as perfect (default) there can only be one module on this layer.
   * Calling AddModuleRow will throw an exception. This is the default behaviour because Gear does
   * not guarantee that the modules are alligned. Displaced modules do not make up a perfect
   * cylinder / plane and have to be treated as separate segments. Finding a neighbouring module/row
   * is not trivial and has to be left to the user or a future Gear version.
   */
		 
  class GearTPCMeasLayer 
    : public TVMeasLayer
  {

  public:
    /** The constructor.
     *  The materials and the type (active or passive) are passed on to the
     *  TVMeasLayer. sigmaX0 [mm] is the constant part of sigmaX, sigmaX1 [mm/sqrt(mm)] 
     *  the z-dependent part, accordingly for sigmaZ.
     *
     *  Module and row have to be specified. They will be added as the first
     *  module/row pair of this measurement layer. 
     *  For a perfect layer modules can be added with AddModuleRow.
     *
     *  Note: This class cannot be instantiated because the parent's geometry dependent
     *  purely virtual
     *  functions like XvToMv are not implemented. This will happen in the cylindrical or planar
     *  implementations.
     *
     *  For inactive layers you will usually leave the sigmas at 0, they have no useful meaning in 
     *  this case.
     */
    GearTPCMeasLayer(TMaterial &min,
		     TMaterial &mout,
		     Int_t      module,
		     Int_t      row,
		     Bool_t     isPerfect,
		     Bool_t     isActive,
		     Double_t   sigmaX0 = 0., //< the constant part of sigmaX
		     Double_t   sigmaX1 = 0., //< the z-dependent part of sigmaX
		     Double_t   sigmaZ0 = 0. , //< the constant part of sigmaZ
		     Double_t   sigmaZ1 = 0.); //< the z-dependent part of sigmaZ
    
    /// The destructor
    virtual  ~GearTPCMeasLayer();

   /**
   * A perfect measurement layer contains all the modules with rows (row segments)
   * that make up the layer.
   */
    virtual std::set< std::pair <int, int> > const & GetModuleRows() const;
    
    /**
     * Add another row on another module which lies on the same cylinder.
     */
    virtual void AddModuleRow(int module, int row);
    
    /**
     * Get the measurement vector (mv) for this layer from a space point (xv)
     */
    virtual TKalMatrix XvToMv    (const TVector3   &xv)   const = 0;

    /**
     * Get the z-depenent resolution in the readout plane 
     * (usually x or r\f$\phi\f$).
     */
    virtual Double_t GetSigmaX(Double_t z) const;

     /**
     * Get the z-depenent resolution in z (drift direction).
     */
   virtual Double_t GetSigmaZ(Double_t z) const;
    
  
    /**
     * Get the flag whether the layer is declared as perfect.
     */
    virtual Bool_t IsPerfect() const;

    /**
     * A virtual function to create the appropriate hit. Depending on the implementation
     * (cylindrical or straight measurement layer) you get the appropriate implementation 
     * of GearTPCHit.
     * It creates a new hit on the heap and hands over the ownership.
     */
    virtual GearTPCHit * createHit(Double_t * meas,
				   Double_t * dmeas,
				   void * hitPointer, 
				   Double_t bField,
				   Double_t vDrift,
				   Int_t           m = kMdim) const = 0;

  protected:
    Double_t fSigmaX0;  // xy resolution
    Double_t fSigmaX1;  // xy resolution
    Double_t fSigmaZ0;  // z  resolution
    Double_t fSigmaZ1;  // z  resolution

    /// A set to hold all the module/row combinations associated to this layer
    std::set< std::pair<int, int> > fModuleRows;

    Bool_t fIsPerfect;
  };

}// namespace kaldet
#endif // GEARTPC_MEASLAYER_H

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("lctpc_gearTPC",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_lctpc_gearTPC_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_lctpc_gearTPC_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_lctpc_gearTPC() {
  TriggerDictionaryInitialization_lctpc_gearTPC_Impl();
}
