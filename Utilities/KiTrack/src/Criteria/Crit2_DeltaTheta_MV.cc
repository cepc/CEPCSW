#include "Criteria/Crit2_DeltaTheta_MV.h"

#include <cmath>
#include <sstream>
//#include <iostream>   // for debugging

using namespace KiTrack;

Crit2_DeltaTheta_MV::Crit2_DeltaTheta_MV ( float deltaThetaMin , float deltaThetaMax ){
   
   
   _deltaThetaMax = deltaThetaMax;  
   _deltaThetaMin = deltaThetaMin;  
   _name = "Crit2_DeltaTheta_MV";
   _type = "2Hit";
   
   _saveValues = false;
   
}



bool Crit2_DeltaTheta_MV::areCompatible( Segment* parent , Segment* child ) {
   
   
   
   if (( parent->getHits().size() == 1 )&&( child->getHits().size() == 1 )){ //a criterion for 1-segments
      

     IHit* a = parent->getHits()[0];
     IHit* b = child->getHits()[0];
     
       
     double thetaa = a->getTheta();
     double thetab = b->getTheta();

     //streamlog_out(DEBUG4) << " theta a = " << (180*thetaa)/M_PI << " theta b = " << (180*thetab)/M_PI << std::endl ;
     
     double deltaTheta = thetaa-thetab;
     
     if (deltaTheta > M_PI) deltaTheta -= 2*M_PI;           //to the range from -pi to pi
     if (deltaTheta < -M_PI) deltaTheta += 2*M_PI;           //to the range from -pi to pi

     float ax = a->getX();
     float ay = a->getY();
     float bx = b->getX();
     float by = b->getY();
     
     if (( by*by + bx*bx < 0.0001 )||( ay*ay + ax*ax < 0.0001 )) deltaTheta = 0.; // In case one of the hits is too close to the origin
      
     deltaTheta = 180.*fabs( deltaTheta ) / M_PI;
     if (_saveValues) _map_name_value["Crit2_DeltaTheta_MV"]= deltaTheta;
     
     
     if ( deltaTheta > _deltaThetaMax ) return false;
     
     if ( deltaTheta < _deltaThetaMin ) return false;
       
     //streamlog_out(DEBUG4) << " delta theta " << deltaTheta << " max " <<  _deltaThetaMax << " min " <<  _deltaThetaMin << std::endl ;    
     
     
   }

   else{
      
      std::stringstream s;
      s << "Crit2_DeltaTheta_MV::This criterion needs 2 segments with 1 hit each, passed was a "
      <<  parent->getHits().size() << " hit segment (parent) and a "
      <<  child->getHits().size() << " hit segment (child).";
      
      
      throw BadSegmentLength( s.str() );
      
      
   }
   
   
   return true;
   
   
   
}


