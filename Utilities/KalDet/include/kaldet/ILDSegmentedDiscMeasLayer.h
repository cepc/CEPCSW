#ifndef __ILDSEGMENTEDDISCMEASLAYER_H__
#define __ILDSEGMENTEDDISCMEASLAYER_H__

/** ILDSegmentedDiscMeasLayer: User defined Segemented Disk Planar KalTest measurement layer class used with ILDPLanarTrackHit. Segments are isosolese trapezoids whose axis of symmetry points to the origin 
 * WARNING: ONLY IMPLEMENTED FOR X AND Y COORDINATES AT FIXED Z
 *
 * @author S.Aplin DESY
 */

#include "TVector3.h"

#include "kaltest/TKalMatrix.h"
#include "kaltest/TPlane.h"
#include "kaltest/KalTrackDim.h"
#include "ILDVMeasLayer.h"

#include "TMath.h"
#include <sstream>
#include <iostream>

#include <vector>

class TVTrackHit;

class ILDSegmentedDiscMeasLayer : public ILDVMeasLayer, public TPlane {
public:
  // Ctors and Dtor
  
  ILDSegmentedDiscMeasLayer(TMaterial &min,
                            TMaterial &mout,
                            double   Bz,
                            double   SortingPolicy,
                            int      nsegments,
                            double   zpos,
                            double   phi0, // defined by the axis of symmerty of the first petal
                            double   trap_rmin,
                            double   trap_height,
                            double   trap_innerBaseLength,
                            double   trap_outerBaseLength,
                            bool     is_active,
                            std::vector<int>      CellIDs,
                            const Char_t    *name = "ILDDiscMeasL");
  
  ILDSegmentedDiscMeasLayer(TMaterial &min,
                            TMaterial &mout,
                            double   Bz,
                            double   SortingPolicy,
                            int      nsegments,
                            double   zpos,
                            double   phi0, // defined by the axis of symmerty of the first petal
                            double   trap_rmin,
                            double   trap_height,
                            double   trap_innerBaseLength,
                            double   trap_outerBaseLength,
                            bool     is_active,
                            const Char_t    *name = "ILDDiscMeasL");
  
  
  // Parrent's pure virtuals that must be implemented
  
  /** Global to Local coordinates */
  virtual TKalMatrix XvToMv    (const TVTrackHit &ht,
                                const TVector3   &xv) const
  { return this->XvToMv(xv); }
  
  
  /** Global to Local coordinates */
  virtual TKalMatrix XvToMv    (const TVector3   &xv) const;
  
  /** Local to Global coordinates */  
  virtual TVector3   HitToXv   (const TVTrackHit &ht) const;
  
  /** Calculate Projector Matrix */
  virtual void       CalcDhDa  (const TVTrackHit &ht,
                                const TVector3   &xv,
                                const TKalMatrix &dxphiada,
                                TKalMatrix &H)  const;
  
  /** Convert LCIO Tracker Hit to an ILDPLanarTrackHit  */
  virtual ILDVTrackHit* ConvertLCIOTrkHit(edm4hep::ConstTrackerHit trkhit) const ;

  /** overloaded version of CalcXingPointWith using closed solution*/
  virtual Int_t    CalcXingPointWith(const TVTrack  &hel,
                                     TVector3 &xx,
                                     Double_t &phi,
                                     Int_t     mode,
                                     Double_t  eps = 1.e-8) const;
  
  /** overloaded version of CalcXingPointWith using closed solution*/
  virtual Int_t    CalcXingPointWith(const TVTrack  &hel,
                                     TVector3 &xx,
                                     Double_t &phi,
                                     Double_t  eps = 1.e-8) const{
    
    return CalcXingPointWith(hel,xx,phi,0,eps);
    
  }
  
  
  
  /** Get the intersection and the CellID, needed for multilayers */
  virtual int getIntersectionAndCellID(const TVTrack  &hel,
                                       TVector3 &xx,
                                       Double_t &phi,
                                       Int_t    &CellID,
                                       Int_t     mode,
                                       Double_t  eps = 1.e-8) const ;

  
  
  /** Check if global point is on surface  */
  virtual Bool_t   IsOnSurface (const TVector3 &xx) const;
  
  /** Get sorting policy for this plane  */
  double GetSortingPolicy() const { return _sortingPolicy; }
  
protected:
  
  double          angular_range_2PI( double phi ) const;
  unsigned int    get_segment_index(double phi) const;
  double          get_segment_phi(unsigned int index) const;
  TVector3        get_segment_centre(unsigned int index) const;
  
private:
  
  double _sortingPolicy;
  int    _nsegments;
  double _trap_rmin;
  double _trap_height;
  double _trap_inner_base_length;
  double _trap_outer_base_length;
  double _trap_tan_beta; // tan of the openning angle of the petal

  double _rmax;
  double _start_phi; // trailing edge of the first sector
  double _segment_dphi;
  
  
  
};



#endif
