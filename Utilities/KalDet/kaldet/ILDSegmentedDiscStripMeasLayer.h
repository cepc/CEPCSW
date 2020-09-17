#ifndef __ILDSEGMENTEDDISCSTRIPMEASLAYER_H__
#define __ILDSEGMENTEDDISCSTRIPMEASLAYER_H__

/** ILDSegmentedDiscStripMeasLayer: User defined Segemented Disk Planar KalTest measurement layer class used with ILDPLanarTrackHit. Segments are isosolese trapezoids whose axis of symmetry points to the origin 
 * WARNING: ONLY IMPLEMENTED FOR X AND Y COORDINATES AT FIXED Z
 *
 * @author S.Aplin DESY
 */
#include "TVector3.h"

#include "kaltest/TKalMatrix.h"
#include "kaltest/TPlane.h"
#include "kaltest/KalTrackDim.h"
#include "ILDSegmentedDiscMeasLayer.h"

#include "TMath.h"
#include <sstream>
#include <iostream>

#include <vector>

class TVTrackHit;

class ILDSegmentedDiscStripMeasLayer : public ILDSegmentedDiscMeasLayer {

public:
  // Ctors and Dtor
  
  ILDSegmentedDiscStripMeasLayer(TMaterial &min,
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
                                 double   stripAngle,
                                 bool     is_active,
                                 std::vector<int>      CellIDs,
                                 const Char_t    *name = "ILDDiscMeasL")
  : ILDSegmentedDiscMeasLayer(min,mout,Bz,SortingPolicy,nsegments,zpos,phi0,trap_rmin,trap_height,trap_innerBaseLength,trap_outerBaseLength,is_active,CellIDs,name), 
  _stripAngle(stripAngle)
  { /* no op */ }
  
  
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
  
private:
  
  double _stripAngle;
  
};



#endif
