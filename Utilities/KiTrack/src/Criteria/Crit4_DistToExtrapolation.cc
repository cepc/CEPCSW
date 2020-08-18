#include "Criteria/Crit4_DistToExtrapolation.h"

#include <cmath>
#include <sstream>

#include "Criteria/SimpleCircle.h"
#include "TVector3.h"



using namespace KiTrack;

Crit4_DistToExtrapolation::Crit4_DistToExtrapolation ( float distMin , float distMax ){
   
   
   _distMin = distMin;
   _distMax = distMax;
   
   _name = "Crit4_DistToExtrapolation";
   _type = "4Hit";
   
   _saveValues = false;
   
}



bool Crit4_DistToExtrapolation::areCompatible( Segment* parent , Segment* child ) {
    
   
   
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
      
         SimpleCircle circle ( ax , ay , bx , by , cx , cy );

         
         double centerX = circle.getCenterX();
         double centerY = circle.getCenterY();
         double R = circle.getRadius();
         
         TVector3 u ( bx - centerX , by - centerY , bz );
         TVector3 v ( cx - centerX , cy - centerY , cz );
         
         
         double deltaPhiParent = v.Phi() - u.Phi(); //angle in xy plane from center of circle, between point 2 and 3
         
         // use this angle and the distance to the next layer to extrapolate
         double zDistParent = fabs( cz - bz );
         double zDistChild  = fabs( dz - cz );
         
         double deltaPhiChild = deltaPhiParent * zDistChild / zDistParent ;
         
         double phiChild = v.Phi() + deltaPhiChild;
         
         double xChildPred = centerX + R* cos(phiChild);
         double yChildPred = centerY + R* sin(phiChild);
         
         
         double DistToPrediction = sqrt ( ( xChildPred- dx )*( xChildPred- dx ) + ( yChildPred- dy )*( yChildPred- dy ) );
         double distNormed = DistToPrediction / zDistChild;   
         
         if (_saveValues) _map_name_value["Crit4_DistToExtrapolation"] = distNormed;
         
         if ( distNormed > _distMax ) return false;
         if ( distNormed < _distMin ) return false;
         
      }
      catch ( InvalidParameter ){
         
         if (_saveValues) _map_name_value["Crit4_DistToExtrapolation"] = -1.;
      
      }
      
   }
   else{
      
      std::stringstream s;
      s << "Crit4_DistToExtrapolation::This criterion needs 2 segments with 3 hits each, passed was a "
      <<  parent->getHits().size() << " hit segment (parent) and a "
      <<  child->getHits().size() << " hit segment (child).";
      
      
      throw BadSegmentLength( s.str() );
      
      
   }
   
   
   return true;
   
   
   
}

