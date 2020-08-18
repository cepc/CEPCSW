#include "Criteria/Crit2_Distance_MV.h"

#include <cmath>
#include <sstream>

using namespace KiTrack;

Crit2_Distance_MV::Crit2_Distance_MV ( float deltaPos2Min , float deltaPos2Max ){
   
   
   _deltaPos2Max = deltaPos2Max;  
   _deltaPos2Min = deltaPos2Min;  
   _name = "Crit2_Distance_MV";
   _type = "2Hit";
   
   _saveValues = false;
   
}



bool Crit2_Distance_MV::areCompatible( Segment* parent , Segment* child ) {
   
   
   
   if (( parent->getHits().size() == 1 )&&( child->getHits().size() == 1 )){ //a criterion for 1-segments
      
      IHit* a = parent->getHits()[0];
      IHit* b = child->getHits()[0];

      float ax = a->getX();
      float ay = a->getY();
      float bx = b->getX();
      float by = b->getY();

      double posa = sqrt(ax*ax + ay*ay);
      double posb = sqrt(bx*bx + by*by);
      
      double deltaPos2 = (posa-posb) * (posa-posb) ;
      
      if (( by*by + bx*bx < 0.0001 )||( ay*ay + ax*ax < 0.0001 )) deltaPos2 = 0.; // In case one of the hits is too close to the origin
      
      if (_saveValues) _map_name_value["Crit2_Distance_MV"]= deltaPos2;
      
        
      if ( deltaPos2 > _deltaPos2Max ) return false;
      
      if ( deltaPos2 < _deltaPos2Min ) return false;
      
      
   }


   else{
      
      std::stringstream s;
      s << "Crit2_Distance_MV::This criterion needs 2 segments with 1 hit each, passed was a "
      <<  parent->getHits().size() << " hit segment (parent) and a "
      <<  child->getHits().size() << " hit segment (child).";
      
      
      throw BadSegmentLength( s.str() );
      
      
   }
   
   
   return true;
   
   
   
}


