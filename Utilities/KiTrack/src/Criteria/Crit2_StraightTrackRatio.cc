#include "Criteria/Crit2_StraightTrackRatio.h"

#include <cmath>
#include <sstream>

using namespace KiTrack;

Crit2_StraightTrackRatio::Crit2_StraightTrackRatio ( float ratioMin, float ratioMax ){
   
   
   _ratioMax = ratioMax; 
   _ratioMin = ratioMin;
   
   _name = "Crit2_StraightTrackRatio";
   _type = "2Hit";
   
   _saveValues = false;
   
}


      
bool Crit2_StraightTrackRatio::areCompatible( Segment* parent , Segment* child ) {
   
   
   
   if (( parent->getHits().size() == 1 )&&( child->getHits().size() == 1 )){ //a criterion for 1-segments
      
      
      
      IHit* a = parent->getHits()[0];
      IHit* b = child-> getHits()[0];
      
      float ax = a->getX();
      float ay = a->getY();
      float az = a->getZ();

      float bx = b->getX();
      float by = b->getY();
      float bz = b->getZ();
      
      //the distance to (0,0) in the xy plane
      double rhoASquared = ax*ax + ay*ay;
      double rhoBSquared = bx*bx + by*by;
      
      

      
      if (_saveValues){
         _map_name_value["Crit2_StraightTrackRatio"]= 1.;
         
      }
     
      
      if( (rhoBSquared >0.) && ( az != 0. ) ){ //prevent division by 0
         
         // the square is used, because it is faster to calculate with the squares than with sqrt, which takes some time!
         double ratioSquared = ( ( rhoASquared * ( bz*bz )  ) / ( rhoBSquared * ( az*az )  ) );
               
         if (_saveValues) _map_name_value["Crit2_StraightTrackRatio"] = sqrt(ratioSquared);
         
         
         if ( ratioSquared > _ratioMax * _ratioMax ) return false;
   
         if ( ratioSquared < _ratioMin * _ratioMin ) return false;
      
      }
      
   }
   else{
      
      std::stringstream s;
      s << "Crit2_StraightTrackRatio::This criterion needs 2 segments with 1 hit each, passed was a "
      <<  parent->getHits().size() << " hit segment (parent) and a "
      <<  child->getHits().size() << " hit segment (child).";
      
      
      throw BadSegmentLength( s.str() );
      
      
   }
   
   
   return true;
   
   
   
}


