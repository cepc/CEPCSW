#ifndef FTDSectorConnector_h
#define FTDSectorConnector_h

#include "KiTrack/ISectorConnector.h"

#include "ILDImpl/SectorSystemFTD.h"



namespace KiTrackMarlin{
   
   /** Used to connect two sectors on the FTD.
    * 
    * 
    * Allows:
    * 
    * - going to layers on the inside (how far see constructor)
    * - going to same petal or petals around (how far see constructor)
    * - jumping to the IP (from where see constructor)
    */   
   class FTDSectorConnector : public ISectorConnector{
      
      
   public:
      
      /** @param layerStepMax the maximum distance of the next layer on the inside
       * 
       *  @param petalStepMax the maximum distance of the next petal (in + and - direction )
       * 
       *  @param lastLayerToIP the highest layer from where a direct jump to the IP is allowed
       */
      FTDSectorConnector ( const SectorSystemFTD* sectorSystemFTD , unsigned layerStepMax , unsigned petalStepMax , unsigned lastLayerToIP);
      
      /** @return a set of all sectors that are connected to the passed sector */
      virtual std::set <int>  getTargetSectors ( int sector );
      
      virtual ~FTDSectorConnector(){};
      
   private:
      
      const SectorSystemFTD* _sectorSystemFTD;
      
      unsigned _layerStepMax;
      unsigned _petalStepMax;
      unsigned _lastLayerToIP;
      
      
   };
   
   
}


#endif

