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
