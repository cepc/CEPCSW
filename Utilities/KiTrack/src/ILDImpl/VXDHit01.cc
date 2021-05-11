#include "ILDImpl/VXDHit01.h"
#include "ILDImpl/SectorSystemVXD.h"

//#include "UTIL/LCTrackerConf.h"
#include <UTIL/ILDConf.h>

#include <iostream>
#include <algorithm>
#include <cmath>
#include <climits>

using namespace KiTrackMarlin;


VXDHit01::VXDHit01( edm4hep::TrackerHit* trackerHit , const SectorSystemVXD* const sectorSystemVXD ){
   
   
   _sectorSystemVXD = sectorSystemVXD;
   
   _trackerHit = trackerHit;

   //Set the position of the VXDHit01
   const edm4hep::Vector3d& pos= trackerHit->getPosition();
   _x = pos[0];
   _y = pos[1]; 
   _z = pos[2]; 


   //find out layer, module, sensor

   UTIL::BitField64  cellID( UTIL::ILDCellID0::encoder_string );

   //cellID.setValue( trackerHit->getCellID0() );
   _layer = cellID[ UTIL::ILDCellID0::layer ] + 1 ;   // + 1 to take into account the IP (considered as layer 0 ) 
   //_layer = cellID[ LCTrackerCellID::layer() ];
   int det_id = 0 ;
   det_id  = cellID[UTIL::ILDCellID0::subdet] ;
   if ( det_id == UTIL::ILDDetID::SIT) { _layer = _layer + 6; }   // need to find a more elegant way...


   double radius = 0;
      
   for (int i=0; i<3; ++i) {
     radius += pos[i]*pos[i];
   }

   radius = sqrt(radius);
      
   double _cosTheta = (pos[2]/radius);
   double _phi = atan2(pos[1],pos[0]);
   double _theta = acos( _cosTheta ) ;
      
   if (_phi < 0.) _phi = _phi + 2*M_PI;   

   // YV, for debugging. Calculate sector here and not through the IVXHit base class
   //calculateSector();

   _sector = _sectorSystemVXD->getSector( _layer, _phi, _cosTheta );

   
   //We assume a real hit. If it is virtual, this has to be set.
   _isVirtual = false;
   
   
}


