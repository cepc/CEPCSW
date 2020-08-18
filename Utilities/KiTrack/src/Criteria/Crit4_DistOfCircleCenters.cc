#include "Criteria/Crit4_DistOfCircleCenters.h"

#include <cmath>
#include <sstream>

#include "Criteria/SimpleCircle.h"


using namespace KiTrack;

Crit4_DistOfCircleCenters::Crit4_DistOfCircleCenters ( float distMin , float distMax ){
   
   
   _distMax = distMax;
   _distMin = distMin;
   
   _name = "Crit4_DistOfCircleCenters";
   _type = "4Hit";
   
   _saveValues = false;
   
   
}



bool Crit4_DistOfCircleCenters::areCompatible( Segment* parent , Segment* child ) {
    
   
   
   if (( parent->getHits().size() == 3 )&&( child->getHits().size() == 3 )){ //this is a criterion for 3-segments
      


      IHit* a = child->getHits()[0];
      IHit* b = child->getHits()[1];
      IHit* c = child->getHits()[2];
      IHit* d = parent-> getHits()[2];

      
      
      float ax = a->getX();
      float ay = a->getY();
//       float az = a->getZ();
      
      float bx = b->getX();
      float by = b->getY();
//       float bz = b->getZ();
      
      float cx = c->getX();
      float cy = c->getY();
//       float cz = c->getZ();
      
      float dx = d->getX();
      float dy = d->getY();
//       float dz = d->getZ();
      
      
      
      try{
      
         SimpleCircle circle1 ( ax , ay , bx , by , cx , cy );
         SimpleCircle circle2 ( bx , by , cx , cy , dx , dy );

         
         float X1 = circle1.getCenterX();
         float Y1 = circle1.getCenterY();

         float X2 = circle2.getCenterX();
         float Y2 = circle2.getCenterY();
         
         
         float distOfCircleCenters = sqrt( (X2-X1)*(X2-X1) + (Y2-Y1)*(Y2-Y1) );
         
         if (_saveValues) _map_name_value["Crit4_DistOfCircleCenters"] = distOfCircleCenters;
         
         
         if ( distOfCircleCenters > _distMax ) return false;
         if ( distOfCircleCenters < _distMin ) return false;
      
      }
      catch ( InvalidParameter ){
         
         if (_saveValues) _map_name_value["Crit4_DistOfCircleCenters"] = 0.;
         
         
      }
      
      
   }
   else{
      
      std::stringstream s;
      s << "Crit4_DistOfCircleCenters::This criterion needs 2 segments with 3 hits each, passed was a "
        <<  parent->getHits().size() << " hit segment (parent) and a "
        <<  child->getHits().size() <<" hit segment (child).";
      
      
      
      throw BadSegmentLength( s.str() );
      
      
   }
   
   return true;
   
   
   
}

