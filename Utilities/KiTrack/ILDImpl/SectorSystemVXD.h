#ifndef SectorSystemVXD_h
#define SectorSystemVXD_h

#include "KiTrack/ISectorSystem.h"

#include <vector>

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
   class SectorSystemVXD : public ISectorSystem{
      
      
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
    SectorSystemVXD( unsigned nLayers , unsigned nDivisionsInPhi , unsigned nDivisionsInTheta );
      

      /** Virtual, because this method is demanded by the Interface ISectorSystem
       * 
       * @return the layer corresponding to the passed sector number
       */
      virtual unsigned getLayer( int sector ) const ;
      
      virtual unsigned getPhi( int sector ) const ;

      virtual unsigned getTheta( int sector ) const ;

      /** @return some information on the sector as string */
      virtual std::string getInfoOnSector( int sector) const;


      /** Calculates the sector number corresponding to the passed parameters
       */
      int getSector( int layer, int phi, int theta ) const ;

      int getSector( int layer, double phi, double cosTheta ) const ;
      
      unsigned getPhiSectors() const ;

      unsigned getThetaSectors() const ;

      unsigned getNLayers() const ;

      virtual ~SectorSystemVXD(){}
      
   private:
      
      int _sectorMax;
      unsigned _nLayers;
      unsigned _nDivisionsInPhi ;
      unsigned _nDivisionsInTheta ;
      
      void checkSectorIsInRange( int sector ) const ;
      
   };



}


#endif


