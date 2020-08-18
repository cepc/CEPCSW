#include "Criteria/Crit3_NoZigZag_MV.h"

#include <cmath>
#include <sstream>

#include "TVector3.h"


using namespace KiTrack;

Crit3_NoZigZag_MV::Crit3_NoZigZag_MV ( float prodMin , float prodMax ){
   
   
   _prodMin = prodMin;
   _prodMax = prodMax;
   
   _name = "Crit3_NoZigZag_MV";
   _type = "3Hit";
   
   _saveValues = false;
   
}



bool Crit3_NoZigZag_MV::areCompatible( Segment* parent , Segment* child ) {
   
   
   
   if (( parent->getHits().size() == 2 )&&( child->getHits().size() == 2 )){ //this is a criterion for 2-segments
      

      IHit* a = child->getHits()[0];
      IHit* b = child->getHits()[1];
      IHit* c = parent-> getHits()[1];

      double thetaa = a->getTheta();
      double thetab = b->getTheta();      
      double thetac = c->getTheta();    

      
      double angleXY1 = thetac - thetab; //the angle between 2-segments in the xy plane
      double angleXY2 = thetab - thetaa;


      if ( a->isVirtual() ) { angleXY2 = 0 ; }
      
      angleXY1 -= 2*M_PI*floor( angleXY1 /2. /M_PI );    //to the range from 0 to 2pi 
      if (angleXY1 > M_PI) angleXY1 -= 2*M_PI;           //to the range from -pi to pi

      angleXY2 -= 2*M_PI*floor( angleXY2 /2. /M_PI );    //to the range from 0 to 2pi 
      if (angleXY2 > M_PI) angleXY2 -= 2*M_PI;           //to the range from -pi to pi
  
      
      // to grad
      angleXY1 *= 180./M_PI;
      angleXY2 *= 180./M_PI;

      float prod = angleXY1 * angleXY2; // if the direction of curvature stays the same, both anlges have the same sign-> and therefore the product is positive


      //streamlog_out(DEBUG4) << " parent layer " << Layerc << "  theta " << (180*thetac)/M_PI << " child first hit layer " << Layera  << " theta = " << (180*thetaa)/M_PI << " child second hit layer " << Layerb  << " theta = " << (180*thetab)/M_PI << " angleXY1 " << angleXY1 <<  " angleXY2 " << angleXY2 <<  " prod " << prod << std::endl ;
      
      

      if (_saveValues) _map_name_value["Crit3_NoZigZag_MV_angle1"] = angleXY1;
      if (_saveValues) _map_name_value["Crit3_NoZigZag_MV_angle2"] = angleXY2;
      if (_saveValues) _map_name_value["Crit3_NoZigZag_MV"] = prod;
      
      if ( prod < _prodMin ) return false;
      if ( prod > _prodMax ) return false;
      
         
      
   }
   else{
      
      std::stringstream s;
      s << "Crit3_NoZigZag_MV::This criterion needs 2 segments with 2 mini-vectors each, passed was a "
      <<  parent->getHits().size() << " hit segment (parent) and a "
      <<  child->getHits().size() << " hit segment (child).";
      
      
      throw BadSegmentLength( s.str() );
      
      
   }
   
   return true;
   
   
   
}
