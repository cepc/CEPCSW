#include "Criteria/Crit3_ChangeRZRatio.h"

#include <cmath>
#include <sstream>

using namespace KiTrack;


Crit3_ChangeRZRatio::Crit3_ChangeRZRatio( float minChange , float maxChange ){
   
   
   _ratioChangeMaxSquared = maxChange*maxChange;
   _ratioChangeMinSquared = minChange*minChange;
   
   _name = "Crit3_ChangeRZRatio";
   _type = "3Hit";
   
   _saveValues = false;
   
}



bool Crit3_ChangeRZRatio::areCompatible( Segment* parent , Segment* child ) {
   
   
   
   if (( parent->getHits().size() == 2 )&&( child->getHits().size() == 2 )){ //a criterion for 2-segments


      IHit* a = child->getHits()[0];
      IHit* b = child->getHits()[1];
      IHit* c = parent-> getHits()[1];
      
      float ax = a->getX();
      float ay = a->getY();
      float az = a->getZ();
      
      float bx = b->getX();
      float by = b->getY();
      float bz = b->getZ();
      
      float cx = c->getX();
      float cy = c->getY();
      float cz = c->getZ();



      // The rz ratios squared
      
      double ratioSquaredParent = 0.; 
      if ( az-bz  != 0. ) ratioSquaredParent = ( (ax-bx)*(ax-bx) + (ay-by)*(ay-by) + (az-bz)*(az-bz) ) / ( (az-bz) * ( az-bz ) );
      
      double ratioSquaredChild = 0.; 
      if ( cz-bz  != 0. ) ratioSquaredChild = ( (cx-bx)*(cx-bx) + (cy-by)*(cy-by) + (cz-bz)*(cz-bz) ) / ( (cz-bz) * ( cz-bz ) );

      double ratioOfRZRatioSquared = 0.;
      
      if (ratioSquaredChild != 0.) ratioOfRZRatioSquared = ratioSquaredParent / ratioSquaredChild;
      
      if (_saveValues) {
         
         _map_name_value["Crit3_ChangeRZRatio_ratioOfRZRatioSquared"] =  ratioOfRZRatioSquared;
         _map_name_value["Crit3_ChangeRZRatio"] = sqrt( ratioOfRZRatioSquared );
         
      }

      if ( ratioOfRZRatioSquared > _ratioChangeMaxSquared ) return false;
      if ( ratioOfRZRatioSquared < _ratioChangeMinSquared ) return false;


   }
   else{
      
      std::stringstream s;
      s << "Crit3_ChangeRZRatio::This criterion needs 2 segments with 2 hits each, passed was a "
      <<  parent->getHits().size() << " hit segment (parent) and a "
      <<  child->getHits().size() << " hit segment (child).";
      
      
      throw BadSegmentLength( s.str() );
      
      
   }
   
   
   
   return true;
   
   
}
