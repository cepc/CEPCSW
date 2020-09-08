#ifndef __ILDSEGMENTEDDISCMEASLAYER_H__
#define __ILDSEGMENTEDDISCMEASLAYER_H__

/** ILDPolygonBarrelMeasLayer: User defined Polygonal Barrel KalTest measurement layer class to be used only for dead material. Segments are planes parallel to the z axis
 *
 *   NOTE: ALL METHODS INVOLVING HITS ARE DISABLED AND CALL EXIT(1)
 *         THIS CLASS IS ONLY MEANT FOR DEAD MATERIAL
 *
 * @author S.Aplin DESY
 */

#include "TVector3.h"

#include "kaltest/TKalMatrix.h"
#include "kaltest/TVSurface.h"
#include "kaltest/KalTrackDim.h"
#include "ILDVMeasLayer.h"

#include "TMath.h"
#include <sstream>

#include "ILDParallelPlanarMeasLayer.h"

#include <vector>

class TVTrackHit;


class ILDPolygonBarrelMeasLayer : public ILDVMeasLayer, public TVSurface {
public:
  // Ctors and Dtor
  

  
  ILDPolygonBarrelMeasLayer(TMaterial &min,
                            TMaterial &mout,
                            double   Bz,
                            double   SortingPolicy,
                            double   r0,       // min distance to the z-axis
                            double   lhalf,    // half length
                            int      nsides,   
                            double   zpos,     // z of the centre 
                            double   phi0,     // phi of the first normal following the xaxis positive rotation
                            std::vector<int>      CellIDs,
                            bool     is_active,
                            const Char_t    *name = "ILDPolygonBarrelMeasL");
  

  ~ILDPolygonBarrelMeasLayer(){ delete _enclosing_cylinder;}
  
  
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

  
  bool IsOutside(const TVector3 &xx) const;
  
  double CalcS(const TVector3 &xx) const;
  
  TMatrixD CalcDSDx(const TVector3 &xx) const;
  
  /** Check if global point is on surface  */
  inline virtual Bool_t   IsOnSurface (const TVector3 &xx) const;
  
  /** Get sorting policy for this plane  */
  double GetSortingPolicy() const { return _sortingPolicy; }
  
private:
  
  double angular_range_2PI( double phi ) const;
  
  unsigned int    get_plane_index(double phi) const;
  
  double _sortingPolicy;

  double   _r0;       // min distance to the z-axis
  double   _lhalf;    // half length
  int      _nsides;   
  double   _zpos;     // z of the centre 
  double   _phi0;     // phi of the first normal following the xaxis positive rotation
  
  double   _segment_dphi;
  double   _start_phi;
  double   _rmax;
  
  std::vector<ILDParallelPlanarMeasLayer> _planes;
  TCylinder* _enclosing_cylinder;
};



#endif
