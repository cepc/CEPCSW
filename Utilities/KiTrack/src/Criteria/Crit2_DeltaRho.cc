#include "Criteria/Crit2_DeltaRho.h"

#include <cmath>
#include <sstream>

using namespace KiTrack;

Crit2_DeltaRho::Crit2_DeltaRho ( float deltaRhoMin , float deltaRhoMax ){
   
   
   _deltaRhoMax = deltaRhoMax;  
   _deltaRhoMin = deltaRhoMin;
   
   _name = "Crit2_DeltaRho";
   _type = "2Hit";
   
   _saveValues = false;
   
}


      
bool Crit2_DeltaRho::areCompatible( Segment* parent , Segment* child ) {
   
   
   
   if (( parent->getHits().size() == 1 )&&( child->getHits().size() == 1 )){ //a criterion for 1-segments
      
      
      
      IHit* a = parent->getHits()[0];
      IHit* b = child-> getHits()[0];
      
      float ax = a->getX();
      float ay = a->getY();

      float bx = b->getX();
      float by = b->getY();
      
      //the distance to (0,0) in the xy plane
      float rhoA =  sqrt( ax*ax + ay*ay );
      float rhoB =  sqrt( bx*bx + by*by );
      
      float deltaRho = rhoA - rhoB;
      
      //first check, if the distance to (0,0) rises --> such a combo could not reach the IP
      if (_saveValues){
         _map_name_value["Crit2_DeltaRho_rhoParent"] = rhoA;
         _map_name_value["Crit2_DeltaRho_rhoChild"] = rhoB;
         _map_name_value["Crit2_DeltaRho"] = deltaRho;
      }
      
      
     
      if ( deltaRho > _deltaRhoMax ) return false;
      if ( deltaRho < _deltaRhoMin ) return false;
      
      
            
   }
   else{
      
      std::stringstream s;
      s << "Crit2_DeltaRho::This criterion needs 2 segments with 1 hit each, passed was a "
      <<  parent->getHits().size() << " hit segment (parent) and a "
      <<  child->getHits().size() << " hit segment (child).";
      
      
      throw BadSegmentLength( s.str() );
      
      
   }
   
   
   return true;
   
   
   
}


