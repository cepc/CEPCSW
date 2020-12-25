#ifndef __ILDVMEASLAYER__
#define __ILDVMEASLAYER__

/** ILDVMeasLayer: Virtual measurement layer class used by ILD[X]MeasLayer Classes.
 *
 * @author S.Aplin DESY
 */

#include "TVector3.h"
#include "kaltest/TKalMatrix.h"
#include "kaltest/TCylinder.h"
#include "kaltest/TVMeasLayer.h"
#include "kaltest/TAttDrawable.h"
#include "kaltest/KalTrackDim.h"
#include "TString.h"
#include "edm4hep/TrackerHitConst.h"

#include <vector>

class TVTrackHit;
class TNode;
class ILDVTrackHit;

namespace edm4hep{
  class TrackerHit;
}

class ILDVMeasLayer : public TVMeasLayer {
public:
  
  static Bool_t kActive;
  static Bool_t kDummy;
  
  /** Get the layer ID */
  inline int getLayerID() const { return _layerID ; } 
  
  /** Get the Cell ID associated with this measurement layer */
  inline const std::vector<int>& getCellIDs() const { return _cellIDs ; }
  
  /** Get the number of Cell ID associated with this measurement layer */
  inline unsigned int getNCellIDs() const { return _cellIDs.size() ; }
    
  /** Get the Magnetic field at the measurement surface */
  inline Double_t GetBz() const { return _Bz; }
  
  /** Convert LCIO Tracker Hit to an ILDPLanarTrackHit  */
  virtual ILDVTrackHit* ConvertLCIOTrkHit(edm4hep::ConstTrackerHit trkhit) const = 0 ;
  
  /** Check whether the measurement layer represents a series of detector elements */
  bool isMultilayer() const { return _isMultiLayer; } 
  
  /** Get the intersection and the CellID, needed for multilayers */
  virtual int getIntersectionAndCellID(const TVTrack  &hel,
                               TVector3 &xx,
                               Double_t &phi,
                               Int_t    &CellID,
                               Int_t     mode,
                               Double_t  eps = 1.e-8) const = 0 ; 
  
protected:
  
  ILDVMeasLayer(TMaterial &min,
                TMaterial &mout,
                Double_t  Bz,
                Bool_t    is_active = ILDVMeasLayer::kActive,
                int CellID = -1 , 
                const Char_t    *name = "ILDMeasL");
  
  ILDVMeasLayer(TMaterial &min,
                TMaterial &mout,
                Double_t  Bz,
                const std::vector<int>& cellIDs,
                Bool_t    is_active = ILDVMeasLayer::kActive,
                const Char_t    *name = "ILDMeasL");
  
  
  
  Double_t _Bz ;       // Magnitude of B-Field in Z
  int _layerID ;
  std::vector<int> _cellIDs ;

  bool _isMultiLayer;
  
private:
  
  
};

#endif
