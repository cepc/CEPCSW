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
