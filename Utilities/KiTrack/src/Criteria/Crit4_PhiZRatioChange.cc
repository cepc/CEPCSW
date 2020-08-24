#include "Criteria/Crit4_PhiZRatioChange.h"

#include <cmath>
#include <sstream>

#include "Criteria/SimpleCircle.h"
#include "TVector3.h"



using namespace KiTrack;

Crit4_PhiZRatioChange::Crit4_PhiZRatioChange ( float changeMin , float changeMax ){
   
   
   _changeMin = changeMin;
   _changeMax = changeMax;
   
   _name = "Crit4_PhiZRatioChange";
   _type = "4Hit";
   
   _saveValues = false;
   
}



bool Crit4_PhiZRatioChange::areCompatible( Segment* parent , Segment* child ) {
    
   
   
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
      float bz = b->getZ();
      
      float cx = c->getX();
      float cy = c->getY();
      float cz = c->getZ();
      
      float dx = d->getX();
      float dy = d->getY();
      float dz = d->getZ();
      
      try{
      
         SimpleCircle circle1 ( ax , ay , bx , by , cx , cy );
         SimpleCircle circle2 ( bx , by , cx , cy , dx , dy );
         
         
         float X1 = circle1.getCenterX();
         float Y1 = circle1.getCenterY();
         float X2 = circle2.getCenterX();
         float Y2 = circle2.getCenterY();
         
         
         
         TVector3 u ( bx - X1, by - Y1 , 0.); //vector from center of circle to point
         TVector3 v ( cx - X1, cy - Y1 , 0.);
         float zDist1 = fabs( cz - bz );
         float phi1 = u.Angle( v );
         float phiZRatio1 = phi1 / zDist1;
         
         TVector3 s ( cx - X2, cy - Y2 , 0.); //vector from center of circle to point
         TVector3 t ( dx - X2, dy - Y2 , 0.);
         float zDist2 = fabs( dz - cz );
         float phi2 = s.Angle( t );
         float phiZRatio2 = phi2 / zDist2;
         
         
         float ratioOfPhiZRatio = phiZRatio1 / phiZRatio2;
         
         if (_saveValues) _map_name_value["Crit4_PhiZRatioChange"] = ratioOfPhiZRatio;
         
         
            
         if ( ratioOfPhiZRatio > _changeMax ) return false;    
         if ( ratioOfPhiZRatio < _changeMin ) return false;
         
      }
      catch ( InvalidParameter ){
       
         if (_saveValues) _map_name_value["Crit4_PhiZRatioChange"] = 1.;
         
      }
      
         
      
   }
   else{
      
      std::stringstream s;
      s << "Crit4_PhiZRatioChange::This criterion needs 2 segments with 3 hits each, passed was a "
      <<  parent->getHits().size() << " hit segment (parent) and a "
      <<  child->getHits().size() << " hit segment (child).";
      
      
      throw BadSegmentLength( s.str() );
      
      
   }
   
   return true;
   
   
   
}

