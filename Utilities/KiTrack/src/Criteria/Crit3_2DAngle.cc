#include "Criteria/Crit3_2DAngle.h"

#include <cmath>
#include <sstream>

using namespace KiTrack;

Crit3_2DAngle::Crit3_2DAngle ( float angleMin, float angleMax ){
   
   
   _cosAngleMin = cos ( angleMax * M_PI / 180. );
   _cosAngleMax = cos ( angleMin * M_PI / 180. );
   
   _name = "Crit3_2DAngle";
   _type = "3Hit";
   
   _saveValues = false;
   
}


      
bool Crit3_2DAngle::areCompatible( Segment* parent , Segment* child ) {
   
   
   //this is not written very beautiful, because this code gets called often and needs to be fast.
   //But it's just a simple angle calculation of two vectors cos(alpha) = u*v/|u||v|
   //
   //Because of speed, I avoided using stuff like sqrt or cos or something in here.
   //That's why it may look a bit odd.
   
   
   
   
   if (( parent->getHits().size() == 2 )&&( child->getHits().size() == 2 )){ //this is a criterion for 2-segments
      


      IHit* a = child->getHits()[0];
      IHit* b = child->getHits()[1];
      IHit* c = parent-> getHits()[1];
      
      float ax = a->getX();
      float ay = a->getY();
      
      float bx = b->getX();
      float by = b->getY();
      
      float cx = c->getX();
      float cy = c->getY();
      
      
      float ux = bx - ax;
      float uy = by - ay;

      
      float vx = cx - bx;
      float vy = cy - by;

   
      //In the numerator there is the vector product of u and v   
      double numerator= ux*vx + uy*vy;
      
      //In the denominator there are the lengths of u and v (here squared)
      double uSquared= ux*ux + uy*uy;
      double vSquared= vx*vx + vy*vy;
      

      
      double denomSquared = uSquared * vSquared;
      
      if (_saveValues){
         
         _map_name_value["Crit3_2DAngle_cos2DAngleSquared"] =  1.;
         _map_name_value["Crit3_2DAngle"] = 0.;
         
      }
      
      if ( denomSquared > 0.){ //don't divide by 0
      
         double cosThetaSquared = numerator * numerator / ( uSquared * vSquared );
         if( cosThetaSquared > 1. ) cosThetaSquared = 1; // prevent rounding errors: cosTheta can mathematically never be bigger than 1!!!
         
         if (_saveValues){
            
            _map_name_value["Crit3_2DAngle_cos2DAngleSquared"] =  cosThetaSquared;
            _map_name_value["Crit3_2DAngle"] = acos( sqrt( cosThetaSquared ) ) * 180. / M_PI;
            
         }
         
         if (cosThetaSquared < _cosAngleMin*_cosAngleMin) return false;
         if (cosThetaSquared > _cosAngleMax*_cosAngleMax) return false;
      
      }
      
      
      

      
      
   }
   else{
      
      std::stringstream s;
      s << "Crit3_2DAngle::This criterion needs 2 segments with 2 hits each, passed was a "
      <<  parent->getHits().size() << " hit segment (parent) and a "
      <<  child->getHits().size() << " hit segment (child).";
      
      
      throw BadSegmentLength( s.str() );
      
      
   }
   
   
   return true;
   
   
   
}

