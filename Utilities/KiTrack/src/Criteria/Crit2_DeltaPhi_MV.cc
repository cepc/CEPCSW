#include "Criteria/Crit2_DeltaPhi_MV.h"

#include <cmath>
#include <sstream>

using namespace KiTrack;

Crit2_DeltaPhi_MV::Crit2_DeltaPhi_MV ( float deltaPhiMin , float deltaPhiMax ){
   
   
   _deltaPhiMax = deltaPhiMax;  
   _deltaPhiMin = deltaPhiMin;  
   _name = "Crit2_DeltaPhi_MV";
   _type = "2Hit";
   
   _saveValues = false;
   
}



bool Crit2_DeltaPhi_MV::areCompatible( Segment* parent , Segment* child ) {
   
   
   
   if (( parent->getHits().size() == 1 )&&( child->getHits().size() == 1 )){ //a criterion for 1-segments
      
      IHit* a = parent->getHits()[0];
      IHit* b = child->getHits()[0];

      double phia = a->getPhi();
      double phib = b->getPhi();
      
      double deltaPhi = phia-phib;
      
      if (deltaPhi > M_PI) deltaPhi -= 2*M_PI;           //to the range from -pi to pi
      if (deltaPhi < -M_PI) deltaPhi += 2*M_PI;           //to the range from -pi to pi
      
      float ax = a->getX();
      float ay = a->getY();
      float bx = b->getX();
      float by = b->getY();
      
      if (( by*by + bx*bx < 0.0001 )||( ay*ay + ax*ax < 0.0001 )) deltaPhi = 0.; // In case one of the hits is too close to the origin
      
      deltaPhi = 180.*fabs( deltaPhi ) / M_PI;
      if (_saveValues) _map_name_value["Crit2_DeltaPhi_MV"]= deltaPhi;
      
        
      if ( deltaPhi > _deltaPhiMax ) return false;
      
      if ( deltaPhi < _deltaPhiMin ) return false;
      
      
   }


   else{
      
      std::stringstream s;
      s << "Crit2_DeltaPhi_MV::This criterion needs 2 segments with 1 hit each, passed was a "
      <<  parent->getHits().size() << " hit segment (parent) and a "
      <<  child->getHits().size() << " hit segment (child).";
      
      
      throw BadSegmentLength( s.str() );
      
      
   }
   
   
   return true;
   
   
   
}


