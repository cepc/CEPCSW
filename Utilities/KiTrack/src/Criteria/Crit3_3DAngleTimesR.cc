#include "Criteria/Crit3_3DAngleTimesR.h"

#include <cmath>
#include <sstream>

#include "Criteria/SimpleCircle.h"


using namespace KiTrack;

Crit3_3DAngleTimesR::Crit3_3DAngleTimesR ( float angleMin, float angleMax ){
   
   
   _angleMin = angleMin;
   _angleMax = angleMax;
   
   _name = "Crit3_3DAngleTimesR";
   _type = "3Hit";
   
   _saveValues = false;
   
}


      
bool Crit3_3DAngleTimesR::areCompatible( Segment* parent , Segment* child ) {
   
   
   
   
   
   
   if (( parent->getHits().size() == 2 )&&( child->getHits().size() == 2 )){ //this is a criterion for 2-segments
      
      
      
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
      
      
      float ux = bx - ax;
      float uy = by - ay;
      float uz = bz - az;
      
      float vx = cx - bx;
      float vy = cy - by;
      float vz = cz - bz;
   
      //In the numerator there is the vector product of u and v   
      double numerator= ux*vx + uy*vy + uz*vz;
      
      //In the denominator there are the lengths of u and v (here squared)
      double uSquared= ux*ux + uy*uy + uz*uz;
      double vSquared= vx*vx + vy*vy + vz*vz;
      

      
      double denomSquared = uSquared * vSquared;
      
      if (_saveValues){
         
         _map_name_value["Crit3_3DAngleTimesR"] = 0.;
         
      }
      
      if ( denomSquared > 0.){ //don't divide by 0
         
         double cosThetaSquared = numerator * numerator / ( uSquared * vSquared );
         if( cosThetaSquared > 1. ) cosThetaSquared = 1; // prevent rounding errors: cosTheta can mathematically never be bigger than 1!!!
         double angle = acos( sqrt( cosThetaSquared ) ) * 180. / M_PI;
         
         
         try{
            
            SimpleCircle circle ( ax , ay , bx , by , cx , cy );
            
            double R = circle.getRadius();
            
            double angleTimesR = R * angle;
            
            if (_saveValues){
               
               _map_name_value["Crit3_3DAngleTimesR"] = angleTimesR;
               
            }
            
            if ( angleTimesR < _angleMin ) return false;
            if ( angleTimesR > _angleMax ) return false;
            
         }
         catch( InvalidParameter ){
            
            
         }
         
         
      }
      
      
      
      
      
   }
   else{
      
      std::stringstream s;
      s << "Crit3_3DAngleTimesR::This criterion needs 2 segments with 2 hits each, passed was a "
      <<  parent->getHits().size() << " hit segment (parent) and a "
      <<  child->getHits().size() << " hit segment (child).";
      
      
      throw BadSegmentLength( s.str() );
      
      
   }
   
   
   return true;
   
   
   
}

