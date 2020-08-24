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
