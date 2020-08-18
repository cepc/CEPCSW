#include "Criteria/Crit3_PT_MV.h"

#include <cmath>
#include <sstream>

#include "Criteria/SimpleCircle.h"

using namespace KiTrack;


Crit3_PT_MV::Crit3_PT_MV( float ptMin , float ptMax , float Bz ){
   
   _ptMin = ptMin;
   _ptMax = ptMax;
   _Bz = Bz;
   
   _name = "Crit3_PT_MV";
   _type = "3Hit";
   
   _saveValues = false;
   
}



bool Crit3_PT_MV::areCompatible( Segment* parent , Segment* child ) {
   
   
   
   if (( parent->getHits().size() == 2 )&&( child->getHits().size() == 2 )){ //a criterion for 2-segments


      IHit* a = child->getHits()[0];
      IHit* b = child->getHits()[1];
      IHit* c = parent-> getHits()[1];
      
      float ax = a->getX();
      float ay = a->getY();
     
      float bx = b->getX();
      float by = b->getY();
      
      float cx = c->getX();
      float cy = c->getY();


      try{
      
         SimpleCircle circle ( ax , ay , bx , by , cx , cy );
      

         double R = circle.getRadius();
         
         
         // check if pt is bigger than _ptMin
         //
         // |omega| = K*Bz/pt
         // R = pt / (K*Bz)
         // pt = R * K *Bz
         //
               
         const double K= 0.00029979; //K depends on the used units
         
         double pt = R * K * _Bz;
            
         if (_saveValues) _map_name_value["Crit3_PT_MV"] =  pt;
               
         if ( pt < _ptMin ) return false;
         if ( pt > _ptMax ) return false;

         
      }
      catch ( InvalidParameter ){
         
         if (_saveValues) _map_name_value["Crit3_PT_MV"] =  0.;
         
      }



   }
   else{
      
      std::stringstream s;
      s << "Crit3_PT_MV::This criterion needs 2 segments with 2 hits each, passed was a "
      <<  parent->getHits().size() << " hit segment (parent) and a "
      <<  child->getHits().size() << " hit segment (child).";
      
      
      throw BadSegmentLength( s.str() );
      
      
   }
   
   
   return true;
   
   
}
