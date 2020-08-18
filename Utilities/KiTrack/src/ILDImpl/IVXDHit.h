#ifndef IVXDHit_h
#define IVXDHit_h

#include <iostream>

#include "edm4hep/TrackerHit.h"
//#include "lcio.h"

#include "KiTrack/IHit.h"

#include "ILDImpl/SectorSystemVXD.h"

namespace KiTrackMarlin{
  /** An interface for a hit for the ILD using an lcio TrackerHit as basis.
   * 
   * It comes along with a layer, phi and theta.
   */   
  class IVXDHit : public IHit{
  public:
    
    edm4hep::TrackerHit* getTrackerHit() { return _trackerHit; };
            
    int getTheta() { return _theta; }
    unsigned getPhi() { return _phi; }
    
    //void setLayer( unsigned layer ){ _layer = layer; calculateSector();}
    //void setPhi( unsigned phi ){ _phi = phi; calculateSector();}
    //void setTheta( unsigned theta ){ _theta = theta; calculateSector();}
    void setLayer( unsigned layer ){ _layer = layer; }
    void setPhi( unsigned phi ){ _phi = phi; }
    void setTheta( unsigned theta ){ _theta = theta; }    
    
    
    virtual const ISectorSystem* getSectorSystem() const { return _sectorSystemVXD; };
    
  protected:
    
    edm4hep::TrackerHit* _trackerHit;
      
    int _layer;
    int _phi;
    int _theta;
    
    const SectorSystemVXD* _sectorSystemVXD;
    
    /** Calculates and sets the sector number
     */
    
    //void calculateSector(){ _sector = _sectorSystemVXD->getSector( _layer, _phi, _theta ); }
  };
}
#endif

