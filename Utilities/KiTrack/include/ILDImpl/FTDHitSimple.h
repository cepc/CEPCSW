#ifndef FTDHitSimple_h
#define FTDHitSimple_h

#include "KiTrack/IHit.h"

#include "ILDImpl/SectorSystemFTD.h"



namespace KiTrackMarlin{
   
   
   /** A hit 
    */   
   class FTDHitSimple : public IHit{
      
      
   public:
      
      FTDHitSimple( float x , float y , float z , int side, unsigned layer , unsigned module, unsigned sensor, const SectorSystemFTD* const sectorSystemFTD );
      
      
      
      virtual const ISectorSystem* getSectorSystem() const { return _sectorSystemFTD; };
      
      virtual ~FTDHitSimple(){}
      
   private:
      
      int _side;
      unsigned _layer;
      unsigned _module;
      unsigned _sensor;
      
      const SectorSystemFTD* _sectorSystemFTD;
      
      void calculateSector(){ _sector = _sectorSystemFTD->getSector( _side, _layer , _module , _sensor ); }
      
   };
   
}


#endif

