#include "ILDImpl/VXDSectorConnector.h"


using namespace KiTrackMarlin;


// Constructor
VXDSectorConnector::VXDSectorConnector( const SectorSystemVXD* sectorSystemVXD , unsigned layerStepMax, unsigned lastLayerToIP, int neighPhi, int neighTheta, int layerMax ){
   
   _sectorSystemVXD = sectorSystemVXD ;
   _layerStepMax = layerStepMax ;
   _lastLayerToIP = lastLayerToIP ;

   _nLayers = sectorSystemVXD->getNLayers();
   _nDivisionsInPhi = sectorSystemVXD->getPhiSectors();
   _nDivisionsInTheta = sectorSystemVXD->getThetaSectors();
   _neighPhi = neighPhi ;
   _neighTheta = neighTheta ;
   _layerMax = layerMax ;
}



std::set< int > VXDSectorConnector::getTargetSectors ( int sector ){
   
   
   std::set <int> targetSectors;

   //std::cout << " check : n layers " << _nLayers << " n phi divisions " << _nDivisionsInPhi << " n theta divisions " << _nDivisionsInTheta << std::endl ;

   // Decode the sector integer,  and take the layer, phi and theta bin
   
   int iTheta = sector/(_nLayers*_nDivisionsInPhi) ;
   
   int iPhi = ((sector - (iTheta*_nLayers*_nDivisionsInPhi)) / _nLayers) ;
   
   int layer = sector - (iTheta*_nLayers*_nDivisionsInPhi) - (iPhi*_nLayers) ; 

   // search for sectors at the neighbouring theta nad phi bins

   int iPhi_Up    = iPhi + _neighPhi;
   int iPhi_Low   = iPhi - _neighPhi;
   int iTheta_Up  = iTheta + _neighTheta; 
   int iTheta_Low = iTheta - _neighTheta;
   if (iTheta_Low < 0) iTheta_Low = 0;
   if (iTheta_Up  >= _nDivisionsInTheta) iTheta_Up = _nDivisionsInTheta-1;
   
   //*************************************************************************************

   //streamlog_out(DEBUG2) << " checking sector " << sector << " layer " << layer << " phi " << iPhi << " theta " << iTheta <<std::endl ;
   
   for( unsigned layerStep = 1; layerStep <= _layerStepMax; layerStep++ ){
     
     if ( layer >= layerStep ){ // +1 makes sense if I use IP as innermost layer
       
       unsigned layerTarget = layer - layerStep;

        if (layerTarget >= 0 && layerTarget < _layerMax ){   // just a test to run cellular automaton over the whole VXD - SIT
	 
	 for (int iPhi = iPhi_Low ; iPhi <= iPhi_Up ; iPhi++){

	   int ip = iPhi;
	   
	   // catch wrap-around
	   if (ip < 0) ip = _nDivisionsInPhi-1;          
	   if (ip >= _nDivisionsInPhi) ip = ip - _nDivisionsInPhi;
	   
	   for (int iTheta = iTheta_Low ; iTheta <= iTheta_Up ; iTheta++){
	     
	     targetSectors.insert( _sectorSystemVXD->getSector ( layerTarget , ip , iTheta ) ); 
	     
	   }
	 }
       }
     }
   }
   

   if ( layer > 0 && ( layer <= _lastLayerToIP ) ){
      
     unsigned layerTarget = 0;
      
     for (int ip = iPhi_Low ; ip <= iPhi_Up ; ip++){
       
       for (int iTheta = iTheta_Low ; iTheta <= iTheta_Up ; iTheta++){
	 
	 //streamlog_out(DEBUG1) << " VXDSectorConnector: from layer " << layer << " to layer " << layerTarget << std::endl ;
	 
	 targetSectors.insert( 0 ) ;
       }
     }
   }
   
										 
   return targetSectors;
   
   
}


