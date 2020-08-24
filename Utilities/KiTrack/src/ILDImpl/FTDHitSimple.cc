#include "ILDImpl/FTDHitSimple.h"


using namespace KiTrackMarlin;

FTDHitSimple::FTDHitSimple( float x , float y , float z , int side, unsigned layer , unsigned module, unsigned sensor, const SectorSystemFTD* const sectorSystemFTD ){
   
   
   _sectorSystemFTD = sectorSystemFTD;
   
   _x = x;
   _y = y; 
   _z = z; 
   
  
   _side   = side;
   _layer  = layer;
   _module = module;
   _sensor = sensor;
   
   
   calculateSector();
   
   
   //We assume a real hit. If it is virtual, this has to be set.
   _isVirtual = false;
   
   
}

