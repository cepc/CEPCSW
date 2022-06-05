#include "ILDImpl/FTDHit01.h"


#include "UTIL/ILDConf.h"

using namespace KiTrackMarlin;


FTDHit01::FTDHit01( edm4hep::TrackerHit trackerHit , const SectorSystemFTD* const sectorSystemFTD ){
   
   
   _sectorSystemFTD = sectorSystemFTD;
   
   _trackerHit = trackerHit;

   //Set the position of the FTDHit01
   const edm4hep::Vector3d& pos= trackerHit.getPosition();
   _x = pos[0];
   _y = pos[1]; 
   _z = pos[2]; 


   //find out layer, module, sensor

   UTIL::BitField64  cellID( UTIL::ILDCellID0::encoder_string );

   cellID.setValue( trackerHit.getCellID() );
     
   _side   = cellID[ UTIL::ILDCellID0::side ];
   _module = cellID[ UTIL::ILDCellID0::module ];
   _sensor = cellID[ UTIL::ILDCellID0::sensor ] - 1;
   _layer = cellID[ UTIL::ILDCellID0::layer ] + 1;
   
   
   calculateSector();
   
   
   //We assume a real hit. If it is virtual, this has to be set.
   _isVirtual = false;
   
   
}


