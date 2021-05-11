#include "ILDImpl/SectorSystemFTD.h"

#include <sstream>

using namespace KiTrackMarlin;


SectorSystemFTD::SectorSystemFTD( unsigned nLayers , unsigned nModules , unsigned nSensors ):
   

_nModules( nModules ),
_nSensors( nSensors ){
   
   _nLayers = nLayers;
   _sectorMax = 2*nLayers*nModules*nSensors - 1;
   
}
   


int SectorSystemFTD::getSector( int side, unsigned layer , unsigned module , unsigned sensor )const {
   
   //check if the values passed are okay:
   if ( ( side!= 1 )&&( side != -1 ) ){
      
      
      std::stringstream s;
      s << "Side has to be either +1 or -1 and not " << side;
      throw OutOfRange( s.str() );
      
   }
   
   if ( layer >= _nLayers ){
      
      std::stringstream s; 
      s << "Layer " << layer << " is too big, the outermost layer is layer " << _nLayers - 1;
      throw OutOfRange( s.str() );
      
   }
   
   if ( module >= _nModules ){
      
      std::stringstream s; 
      s << "Module " << module << " is too big, the highest module is module " << _nModules - 1;
      throw OutOfRange( s.str() );
      
   }
   
   if ( sensor >= _nSensors ){
      
      std::stringstream s;
      s << "Sensor " << sensor << " is too big, the highest sensor is sensor " << _nSensors - 1;
      throw OutOfRange( s.str() );
      
   }   
   
   unsigned multiplicator=1;
   
   
   int sector = sensor;
   multiplicator *= _nSensors; //there are nSensors possible values for sensor
   
   sector += module * multiplicator;
   multiplicator *= _nModules;
   
   sector += layer * multiplicator;
   multiplicator *= _nLayers;
   
   
   sector += ( (side + 1 )/2 ) * multiplicator;             // (side+1) /2 gives 0 for backward (-1) and 1 for forward (+1)
   /*
   streamlog_out( DEBUG0 ) << " Sector of side " << side
   << ", layer " << layer
   << ", module " << module
   << ", sensor " << sensor
   << " == " << sector << "\n";
   */
   
   return sector;
   
}






int SectorSystemFTD::getSide( int sector ) const {
   
   checkSectorIsInRange( sector );
   
   
   
   int side = ( sector / ( _nSensors * _nModules * _nLayers ) ) % 2; //this is an integerdivision --> we will get the floor authomatically
   
   side = side*2 - 1 ; //from 0 and 1 to  -1 and 1
   
//    streamlog_out( DEBUG0 ) << "\n Sector " << sector << " == Side " << side;
   
   return side;
   
}

unsigned SectorSystemFTD::getLayer( int sector ) const {
   
   checkSectorIsInRange( sector );
   
   unsigned layer = ( sector / (  _nSensors * _nModules ) ) % _nLayers; //this is an integerdivision --> we will get the floor authomatically
   
//    streamlog_out( DEBUG0 ) << "\n Sector " << sector << " == Layer " << layer;
   
   return layer;
   
   
}

unsigned SectorSystemFTD::getModule( int sector ) const {
   
   
   checkSectorIsInRange( sector );
   
   unsigned module = ( sector / ( _nSensors ) ) % _nModules; //this is an integerdivision --> we will get the floor authomatically
   
//    streamlog_out( DEBUG0 ) << "\n Sector " << sector << " == Module " << module;
   
   return module;
   
   
   
}
unsigned SectorSystemFTD::getSensor( int sector ) const {
   
   
   checkSectorIsInRange( sector );
   
   unsigned sensor = ( sector ) % _nSensors; 
   
//    streamlog_out( DEBUG0 ) << "\n Sector " << sector << " == Sensor " << sensor;
   
   return sensor;
   
}


void SectorSystemFTD::checkSectorIsInRange( int sector ) const {


   if ( sector > _sectorMax ){
      
      std::stringstream s;
      s << "SectorSystemFTD:\n Sector " 
        << sector << " is too big, the highest possible number for a sector in this configuration of FTDSegRepresentation is"
        << _sectorMax 
        << ".\nThe configuration is: nLayers = " << _nLayers
        << ", nModules = " << _nModules
        << ", nSensors = " << _nSensors 
        << "\n With 2 sides (forward and backward) this gives sectors from 0 to 2*" 
        << _nLayers << "*"  
        << _nModules << "*"
        << _nSensors << " -1 = " << 2*_nLayers*_nModules*_nSensors -1 ;
      throw OutOfRange( s.str() );
      
   }  

}

std::string SectorSystemFTD::getInfoOnSector( int sector ) const{
   
   
   std::stringstream s;
   s << " (si" << getSide(sector)  
     << ",la" << getLayer(sector)
     << ",mo" << getModule(sector) 
     << ",se" << getSensor(sector) 
     << ")";
   
   
   return s.str();   
   
   
}


