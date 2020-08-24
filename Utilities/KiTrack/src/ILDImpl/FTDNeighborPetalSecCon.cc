#include "ILDImpl/FTDNeighborPetalSecCon.h"


using namespace KiTrackMarlin;


FTDNeighborPetalSecCon::FTDNeighborPetalSecCon( const SectorSystemFTD* sectorSystemFTD ){
   
   _sectorSystemFTD = sectorSystemFTD;
     
}



std::set< int > FTDNeighborPetalSecCon::getTargetSectors ( int sector ){
   
   
   
   std::set <int> targetSectors;
   
   
   int side = _sectorSystemFTD->getSide( sector );
   unsigned layer = _sectorSystemFTD->getLayer( sector );
   unsigned petal = _sectorSystemFTD->getModule( sector );
//    unsigned sensor = _sectorSystemFTD->getSensor( sector );
   
//    unsigned nLayers = _sectorSystemFTD->getNumberOfLayers();
   unsigned nPetals = _sectorSystemFTD->getNumberOfModules();
   unsigned nSensors = _sectorSystemFTD->getNumberOfSensors();
   
   
   unsigned petalToTheLeft = petal - 1; //the names left and right are arbitrary, as it of course depends on from where one looks.
   unsigned petalToTheRight = petal + 1;
   
   //Now we have to make sure that we didn't leave the petal range 0 to nPetals-1
   if (petal == 0) petalToTheLeft = nPetals - 1;
   if (petal == nPetals - 1) petalToTheRight = 0;
   
   
   
   for ( unsigned iSensor=0; iSensor < nSensors ; iSensor++ ){ //over all sensors
      
      
      targetSectors.insert( _sectorSystemFTD->getSector ( side , layer , petalToTheLeft , iSensor ) ); 
      targetSectors.insert( _sectorSystemFTD->getSector ( side , layer , petalToTheRight , iSensor ) );
      
   }
   
   
   
   return targetSectors;
   
   
}


