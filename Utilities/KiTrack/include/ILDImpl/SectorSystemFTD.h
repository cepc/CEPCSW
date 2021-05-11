#ifndef SectorSystemFTD_h
#define SectorSystemFTD_h

#include "KiTrack/ISectorSystem.h"

using namespace KiTrack;

namespace KiTrackMarlin{

   /** A Sector System class for the Forward Tracking Disks FTD in the ILD.
    * 
    * It calculates sectors from the side, layer, sensor and module and vice versa.
    * 
    * @param side: +1 for forward, -1 for backward
    * 
    * @param layer: layer of FTD: 0 is the layer of the IP, 1 is the first FTD disk and so on.
    * 
    * @param module: module
    * 
    * @param sensor: the sensor on the module
    * 
    * 
    */ 
   class SectorSystemFTD : public ISectorSystem{
      
      
   public:
      
      /**Constructor
       * 
       * @param nLayers the number of possible layers. The layers from 0 to n-1 will be available. Keep in mind,
       * that layer 0 is used for the IP.
       * 
       * @param nModules the number of modules per disk.
       * 
       * @param nSensors the number of sensors on one module.
       */
      SectorSystemFTD( unsigned nLayers , unsigned nModules , unsigned nSensors );
      
      
      /** Calculates the sector number corresponding to the passed parameters
       */
      int getSector( int side, unsigned layer , unsigned module , unsigned sensor ) const ;
      
      
      /** Virtual, because this method is demanded by the Interface ISectorSystem
       * 
       * @return the layer corresponding to the passed sector number
       */
      virtual unsigned getLayer( int sector ) const ;
      
      
      /** @return some information on the sector as string */
      virtual std::string getInfoOnSector( int sector ) const;
      
      /** @return the side the sector is on (+1 = forward, -1 = backward)
       */
      int getSide( int sector ) const ;
      
      /** @return the module of the sector
       */
      unsigned getModule( int sector ) const ; 
      
      /** @return the sensor of the sector
       */
      unsigned getSensor( int sector ) const ;
      
      
      
      unsigned getNumberOfModules() const { return _nModules; }
      unsigned getNumberOfSensors() const { return _nSensors; }
      
      virtual ~SectorSystemFTD(){}
      
   private:
      
      unsigned _nModules;
      unsigned _nSensors;
      
      int _sectorMax;
      
      void checkSectorIsInRange( int sector ) const ;
      
   };



}


#endif


