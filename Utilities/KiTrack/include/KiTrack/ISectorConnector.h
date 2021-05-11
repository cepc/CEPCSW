#ifndef ISectorConnector_h
#define ISectorConnector_h

#include <set>


namespace KiTrack{
   
   
   /** Abstract Base Class for SectorConnectors.
    * 
    * A SectorConnector is pretty simple: you put in a sector (int) and get back a bunch of other sectors in a set.
    * What it can be used for: Suppose you want to search for the path of a particle through a detector
    * and you know what possible ways it can go, than you can create a SectorConnector to emulate that.
    * So if you have a hit in sector 43 (whatever 43 means) and you think, that from there the particle will
    * only go to either sector 46, 47 or 48, then you would write a SectorConnector that will return
    * a set of 46,47,48 when the method getTargetSectors( 43 ) is called.
    */
   class ISectorConnector{
      
      
   public:
      
      /** @return a set of sectors somehow linked to the passed one */
      virtual std::set <int> getTargetSectors ( int ) = 0;
      
      virtual ~ISectorConnector(){};
      
   };
   
   
   
}



#endif



