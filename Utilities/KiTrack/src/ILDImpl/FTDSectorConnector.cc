#include "ILDImpl/FTDSectorConnector.h"


using namespace KiTrackMarlin;



FTDSectorConnector::FTDSectorConnector( const SectorSystemFTD* sectorSystemFTD , unsigned layerStepMax , unsigned petalStepMax, unsigned lastLayerToIP){
   
   _sectorSystemFTD = sectorSystemFTD;
   _layerStepMax = layerStepMax;
   _lastLayerToIP = lastLayerToIP;
   _petalStepMax = petalStepMax;
   
}



std::set< int > FTDSectorConnector::getTargetSectors ( int sector ){
   
   
   
   std::set <int> targetSectors;
   
   
   int side = _sectorSystemFTD->getSide( sector );
   unsigned layer = _sectorSystemFTD->getLayer( sector );
   unsigned module = _sectorSystemFTD->getModule( sector );
//    unsigned sensor = _sectorSystemFTD->getSensor( sector );
   
//    unsigned nLayers = _sectorSystemFTD->getNumberOfLayers();
   unsigned nModules = _sectorSystemFTD->getNumberOfModules();
   unsigned nSensors = _sectorSystemFTD->getNumberOfSensors();
   
   
   for( unsigned layerStep = 1; layerStep <= _layerStepMax; layerStep++ ){
      
      
      
      if ( layer >= layerStep +1 ){ //other wise the we could jump past layer 1, ( layer 0 is covered below)
         
         
         unsigned layerTarget = layer - layerStep;
         
         
         for ( unsigned iSensor=0; iSensor < nSensors ; iSensor++){ //over all sensors
            
            
            for ( int iPetal= int(module) - _petalStepMax; iPetal <= int(module) + int(_petalStepMax) ; iPetal++ ){ 
               
               //if iPetal is out of the range from 0 to nModules-1, move it back there. 
               //And of course use a different variable for that. 
               //(Or else we would create and endless loop: imagine we have iPetal = 16 and set it back to 0--> the loop will continue from there until it reaches 16 again and so on...)
               int iModule = iPetal;
               while( iModule < 0 ) iModule+= nModules;
               while( iModule >= int(nModules) ) iModule -= nModules;
               
               targetSectors.insert( _sectorSystemFTD->getSector ( side , layerTarget , iModule , iSensor ) ); 
               
            }
            
         }
      
      }
      
   }
   
   //Allow jumping to layer 0 from layer _lastLayerToIP or less
   if ( ( layer >= 1 )&& ( layer <= _lastLayerToIP ) ){
      
      
      unsigned layerTarget = 0;
      
      for ( unsigned iModule=0; iModule < nModules ; iModule++){ //over all modules
         
         for ( unsigned iSensor=0; iSensor < nSensors ; iSensor++ ){ //over all sensors
            
            
            targetSectors.insert( _sectorSystemFTD->getSector ( side , layerTarget , iModule , iSensor ) ); 
            
         }
         
      }
      
   }
   
   return targetSectors;
   
   
}


