#ifndef VXDSectorConnector_h
#define VXDSectorConnector_h

#include "KiTrack/ISectorConnector.h"

#include "ILDImpl/SectorSystemVXD.h"



namespace KiTrackMarlin{
   
   /** Used to connect two sectors on the VXD.
    * 
    * 
    * Allows:
    * 
    * - going to layers on the inside (how far see constructor)
    * - jumping to the IP (from where see constructor)
    */   
   class VXDSectorConnector : public ISectorConnector{
      
      
   public:
      
    VXDSectorConnector ( const SectorSystemVXD* sectorSystemVXD , unsigned layerStepMax, unsigned lastLayerToIP, int neighPhi = 8,  int neighTheta = 1, int layerMax = 10  ) ;
      
      /** @return a set of all sectors that are connected to the passed sector */
      virtual std::set <int>  getTargetSectors ( int sector );
      
      virtual ~VXDSectorConnector(){};
      
   private:
      
      const SectorSystemVXD* _sectorSystemVXD;
      
      unsigned _layerStepMax;
      unsigned _nLayers;
      unsigned _lastLayerToIP;
      unsigned _nDivisionsInPhi ;
      unsigned _nDivisionsInTheta ; 
      int _layerMax ;   
      int _neighTheta ;
      int _neighPhi ;
   };
   
   
}


#endif

