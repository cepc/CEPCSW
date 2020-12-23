#ifndef FTDNeighborPetalSecCon_h
#define FTDNeighborPetalSecCon_h

#include "KiTrack/ISectorConnector.h"

#include "ILDImpl/SectorSystemFTD.h"



namespace KiTrackMarlin{
   
   /** Used to connect two sectors.
    * 
    * Allows:
    * 
    * - Connections to the neighbouring petals (the one to the left and the one to the right on the same layer and side)
    * 
    */   
   class FTDNeighborPetalSecCon : public ISectorConnector{
      
      
   public:
      
      /**
       * 
       */
      FTDNeighborPetalSecCon ( const SectorSystemFTD* sectorSystemFTD );
      
      virtual std::set <int>  getTargetSectors ( int sector );
      
      virtual ~FTDNeighborPetalSecCon(){};
      
   private:
      
      const SectorSystemFTD* _sectorSystemFTD;
      
      
      
   };
   
   
}


#endif

