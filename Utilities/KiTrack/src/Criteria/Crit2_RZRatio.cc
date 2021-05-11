#include "Criteria/Crit2_RZRatio.h"

#include <cmath>
#include <sstream>


using namespace KiTrack;

Crit2_RZRatio::Crit2_RZRatio ( float ratioMin, float ratioMax ){
   
   
   _ratioMax = ratioMax;
   _ratioMin = ratioMin;
   
   _name = "Crit2_RZRatio";
   _type = "2Hit";
   
   _saveValues = false;
   
   
   
}

  
bool Crit2_RZRatio::areCompatible( Segment* parent , Segment* child ) {
   
   
   
   if (( parent->getHits().size() == 1 )&&( child->getHits().size() == 1 )){ //a criterion for 1-segments
      
      
      IHit* a = parent->getHits()[0];
      IHit* b = child-> getHits()[0];
      
      float ax = a->getX();
      float ay = a->getY();
      float az = a->getZ();
      
      float bx = b->getX();
      float by = b->getY();
      float bz = b->getZ();
      
      // the square is used, because it is faster to calculate with the squares than with sqrt, which takes some time!
      double ratioSquared = 0.; 
      if ( az-bz  != 0. ) ratioSquared = ( (ax-bx)*(ax-bx) + (ay-by)*(ay-by) + (az-bz)*(az-bz) ) / ( (az-bz) * ( az-bz ) );
      
      
      if (_saveValues) _map_name_value[ "Crit2_RZRatio"] = sqrt( ratioSquared );
      

      
      if ( ratioSquared > _ratioMax * _ratioMax ) return false;
      if ( ratioSquared < _ratioMin * _ratioMin ) return false;
  
      
      
   }
   else{
      
      std::stringstream s;
      s << "Crit2_RZRatio::This criterion needs 2 segments with 1 hit each, passed was a "
      <<  parent->getHits().size() << " hit segment (parent) and a "
      <<  child->getHits().size() << " hit segment (child).";
      
      
      throw BadSegmentLength( s.str() );
      
      
   }
   
   
   return true;
   
   
   
}






