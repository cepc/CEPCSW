#ifndef IMiniVector_h
#define IMiniVector_h

#include <iostream>

#include "edm4hep/TrackerHit.h"

#include "KiTrack/IHit.h"
#include "ILDImpl/MiniVector.h"

#include "ILDImpl/SectorSystemVXD.h"

namespace KiTrackMarlin{
   
   
   /** An interface for a mini-vector for the ILD 
    */   
   class IMiniVector : public IHit{
      
      
   public:
      
      
      MiniVector* getMiniVector() { return _miniVector; };
      
      
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
      
      MiniVector* _miniVector;
      
      
      int _layer;
      int _phi;
      int _theta;
      
      const SectorSystemVXD* _sectorSystemVXD;
      
      /** Calculates and sets the sector number
       */

      //void calculateSector(){ _sector = _sectorSystemMV->getSector( _layer, _phi, _theta ); }
      
   };
   
}


#endif

