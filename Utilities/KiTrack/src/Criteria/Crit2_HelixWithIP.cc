#include "Criteria/Crit2_HelixWithIP.h"

#include <cmath>
#include <sstream>

#include "Criteria/SimpleCircle.h"

using namespace KiTrack;

Crit2_HelixWithIP::Crit2_HelixWithIP ( float ratioMin , float ratioMax ){
   
   
   _ratioMax = ratioMax;  
   _ratioMin = ratioMin;  
   
   _name = "Crit2_HelixWithIP";
   _type = "2Hit";
   
   _saveValues = false;
   
}


      
bool Crit2_HelixWithIP::areCompatible( Segment* parent , Segment* child ) {
   
   
   
   if (( parent->getHits().size() == 1 )&&( child->getHits().size() == 1 )){ //a criterion for 1-segments
      
      
      
      IHit* a = parent->getHits()[0];
      IHit* b = child-> getHits()[0];
      
      float ax = a->getX();
      float ay = a->getY();
      float az = a->getZ();

      float bx = b->getX();
      float by = b->getY();
      float bz = b->getZ();
      
      
      
      float ratio = 1.;
      
      try{
         
         SimpleCircle circle ( 0. , 0. , ax , ay , bx , by );
         
         float x = circle.getCenterX();
         float y = circle.getCenterY();
         
         
         //calculate the angle in the circle from the IP to the child hit:
         
         //vector from the center of circle to the IP
         float ux = 0. - x;
         float uy = 0. - y;
         
         //vector from the center of circle to the child
         float vx = bx - x;
         float vy = by - y;

         //vector from the center of circle to the parent
         float wx = ax - x;
         float wy = ay - y;
         
         //the angle between u and v (angle between two vectors: cos(alpha) = u*v/|u||v| )
         float alpha1 = acos ( (ux*vx + uy*vy ) / sqrt(( ux*ux + uy*uy ) * ( vx*vx + vy*vy )) );
         
         //the angle between v and w
         float alpha2 = acos ( (vx*wx + vy*wy ) / sqrt(( vx*vx + vy*vy ) * ( wx*wx + wy*wy )) );
         
         float deltaz1 = bz - 0.;
         float deltaz2 = az - bz;
         

         
         if ( ( alpha2 != 0) && (deltaz1 != 0) ) ratio = ( alpha1 * deltaz2 ) / ( alpha2 * deltaz1 );
         


         
      }
      catch ( InvalidParameter ){
         
        
      }
      


      
      if (_saveValues) _map_name_value["Crit2_HelixWithIP"]= ratio;
       
         
      if ( ratio > _ratioMax ) return false;

      if ( ratio < _ratioMin ) return false;
      
      
      
   }
   else{
      
      std::stringstream s;
      s << "Crit2_HelixWithIP::This criterion needs 2 segments with 1 hit each, passed was a "
      <<  parent->getHits().size() << " hit segment (parent) and a "
      <<  child->getHits().size() << " hit segment (child).";
      
      
      throw BadSegmentLength( s.str() );
      
      
   }
   
   
   return true;
   
   
   
}


