#ifndef ISectorSystem_h
#define ISectorSystem_h

#include <string>

#include "KiTrack/KiTrackExceptions.h"

namespace KiTrack{
   
   
   /** An interface for Sector Systems. 
    * 
    * A sector is a code for a place. So it can for example equal a sensor somewhere in a detector or be something
    * abstract like the IP.
    * 
    * A sector system is able to take a sector (integer number) and give back information
    * about it. This can be things like the number of the sensor or the rough distance from the IP or such
    * things. Or even neighbouring sectors. 
    * But this is all dependent on the circumstances of the detectors and their representation. 
    * 
    * What all SectorSystems have in common is the layer: SectorSystems must be able to
    * return the layer of a sector and how many layers there are all in all. 
    */  
   class ISectorSystem{
      
      
   public:
      
      /** @return the layer of the corresponding sector. */
      virtual unsigned getLayer( int sector ) const  =0;
      
      /** @return the number of layers in the sector system. */
      unsigned getNumberOfLayers() const { return _nLayers; };
      
      /** @return some information on the sector as string */
      virtual std::string getInfoOnSector( int sector ) const = 0;
      
      virtual ~ISectorSystem(){}
      
      
      
   protected:
      
      
      unsigned _nLayers{};
      
      
      
   };
   
   
}


#endif

