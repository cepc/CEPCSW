#ifndef IHit_h
#define IHit_h

#include "KiTrack/ISectorSystem.h"

namespace KiTrack{
   
   
   /** An Interface for hits.
    * 
    * They have of course a position, accessible with getX(), getY() and getZ().
    * 
    * Also a hit has a sector. A sector is just an integer, telling where the hit lies. For example
    * a sector 141 could mean 1st layer, 4th module, 1st sensor. But how information is encoded is dealt
    * with in the SectorSystem classes (see ISectorSystem for the Abstract Base Class ).
    */   
   class IHit{
      
      
   public:

       /** @return the x position */
      float getX() const { return _x; }
      
      /** @return the y position */
      float getY() const { return _y; }
      
      /** @return the z position */
      float getZ() const { return _z; }
      
      /** @return the sector, where the hit lies n */
      int getSector() const { return _sector; }
      
      /** @return the used SectorSystem that encodes the sector of the hit*/
      virtual const ISectorSystem* getSectorSystem() const =0;
      
      /** @return the layer of the hit */
      unsigned getLayer() const {return getSectorSystem()->getLayer( _sector ); }
      
      /** @return whether the hit is virtual (not a real hit). */ 
      bool isVirtual() const { return _isVirtual; };
      void setIsVirtual( bool isVirtual ){ _isVirtual = isVirtual; };
      
      /** @return the distance to an other hit */
      float distTo( IHit* otherHit );
      
      /** @return a string containing the position of the hit */
      std::string getPositionInfo();

      /** to be used with mini-vectors: @return 3D angle */  //YV
      double get3DAngle( IHit* /*otherHit*/) const { return _3DAngle ; }

      /** to be used with mini-vectors: @return azimuth angle */  //YV
      double getPhi() const { return _phiMV ; } 

      /** to be used with mini-vectors: @return polar angle */  //YV
      double getTheta() const { return _thetaMV ; } 

      virtual ~IHit(){}
      
   protected:
      
      
      float _x{};
      float _y{};
      float _z{};

      double _3DAngle{};  //YV
      double _phiMV{};     //YV
      double _thetaMV{};     //YV
      
      int _sector{};
      
      
      bool _isVirtual{};
      
      
   };
   
   
   
}


#endif




