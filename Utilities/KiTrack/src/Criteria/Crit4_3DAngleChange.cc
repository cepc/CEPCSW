#include "Criteria/Crit4_3DAngleChange.h"

#include <cmath>
#include <sstream>

#include "TVector3.h"


using namespace KiTrack;

Crit4_3DAngleChange::Crit4_3DAngleChange ( float changeMin , float changeMax ){
   
   
   _changeMax = changeMax;
   _changeMin = changeMin;
   
   _name = "Crit4_3DAngleChange";
   _type = "4Hit";
   
   _saveValues = false;
   
}



bool Crit4_3DAngleChange::areCompatible( Segment* parent , Segment* child ) {
    
      
   if (( parent->getHits().size() == 3 )&&( child->getHits().size() == 3 )){ //this is a criterion for 3-segments
      


      IHit* a = child->getHits()[0];
      IHit* b = child->getHits()[1];
      IHit* c = child->getHits()[2];
      IHit* d = parent-> getHits()[2];
      
      float ax = a->getX();
      float ay = a->getY();
      float az = a->getZ();
      
      float bx = b->getX();
      float by = b->getY();
      float bz = b->getZ();
      
      float cx = c->getX();
      float cy = c->getY();
      float cz = c->getZ();
      
      float dx = d->getX();
      float dy = d->getY();
      float dz = d->getZ();
      

      
      TVector3 outerVec  (bx-ax , by-ay , bz-az );
      TVector3 middleVec (cx-bx , cy-by , cz-bz );
      TVector3 innerVec  (dx-cx , dy-cy , dz-cz );
      
  
      
      
      double angleXY1 = outerVec.Angle( middleVec ); 
      double angleXY2 = middleVec.Angle( innerVec );
      
      angleXY1 -= 2*M_PI*floor( angleXY1 /2. /M_PI );    //to the range from 0 to 2pi 
      if (angleXY1 > M_PI) angleXY1 -= 2*M_PI;           //to the range from -pi to pi

      angleXY2 -= 2*M_PI*floor( angleXY2 /2. /M_PI );    //to the range from 0 to 2pi 
      if (angleXY2 > M_PI) angleXY2 -= 2*M_PI;           //to the range from -pi to pi
  
      
      float ratioOf3DAngles = angleXY1 / angleXY2 ;

      if (_saveValues) _map_name_value["Crit4_3DAngleChange"] = ratioOf3DAngles;
      
      if ( ratioOf3DAngles > _changeMax ) return false;    
      if ( ratioOf3DAngles < _changeMin ) return false;
      
      
      
   }
   else{
      
      std::stringstream s;
      s << "Crit4_3DAngleChange::This criterion needs 2 segments with 3 hits each, passed was a "
      <<  parent->getHits().size() << " hit segment (parent) and a "
      <<  child->getHits().size() << " hit segment (child).";
      
      
      throw BadSegmentLength( s.str() );
      
      
   }
   
   return true;
   
   
   
}

