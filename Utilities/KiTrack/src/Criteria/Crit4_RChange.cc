#include "Criteria/Crit4_RChange.h"

#include <cmath>
#include <sstream>

#include "Criteria/SimpleCircle.h"



using namespace KiTrack;

Crit4_RChange::Crit4_RChange ( float changeMin , float changeMax ){
   
   
   _changeMin = changeMin;
   _changeMax = changeMax;
   
   _name = "Crit4_RChange";
   _type = "4Hit";
   
   _saveValues = false;
   
}



bool Crit4_RChange::areCompatible( Segment* parent , Segment* child ) {
    
   
   
   if (( parent->getHits().size() == 3 )&&( child->getHits().size() == 3 )){ //this is a criterion for 3-segments
      


      IHit* a = child->getHits()[0];
      IHit* b = child->getHits()[1];
      IHit* c = child->getHits()[2];
      IHit* d = parent-> getHits()[2];
      
      float ax = a->getX();
      float ay = a->getY();
      
      float bx = b->getX();
      float by = b->getY();
       
      float cx = c->getX();
      float cy = c->getY();
      
      float dx = d->getX();
      float dy = d->getY();
      
      try{
      
         SimpleCircle circle1 ( ax , ay , bx , by , cx , cy );
         SimpleCircle circle2 ( bx , by , cx , cy , dx , dy );
         
         float R1 = circle1.getRadius();
         float R2 = circle2.getRadius();
         
         float ratioOfR = 1.;
         if (R2 > 0) ratioOfR = R1/R2;
         
         if (_saveValues) _map_name_value["Crit4_RChange"] = ratioOfR;
         
         
            
         if ( ratioOfR > _changeMax ) return false;    
         if ( ratioOfR < _changeMin ) return false;
         
      }
      catch ( InvalidParameter ){
         
         if (_saveValues) _map_name_value["Crit4_RChange"] = 1.;
         
      }
      
   }
   else{
      
      std::stringstream s;
      s << "Crit4_RChange::This criterion needs 2 segments with 3 hits each, passed was a "
      <<  parent->getHits().size() << " hit segment (parent) and a "
      <<  child->getHits().size() << " hit segment (child).";
      
      
      throw BadSegmentLength( s.str() );
      
      
   }
   
   
   return true;
   
   
   
}

