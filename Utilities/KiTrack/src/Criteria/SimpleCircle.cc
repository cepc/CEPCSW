#include "Criteria/SimpleCircle.h"

#include <cmath>
#include <sstream>



using namespace KiTrack;

SimpleCircle::SimpleCircle( double x1 , double y1 , double x2 , double y2 , double x3, double y3 ) {
  
  
   
  // 1. check if they are not in a line, i.e. the slopes are parallel (or two or more points are identical)
  
  if ( (x2 -x1)*(y3 - y2) == (x3 - x2)*(y2 - y1) ){
     
     
     std::stringstream s;
     
     s << "SimpleCircle::The 3 points are on one line in xy-space: x1 = "
       <<  x1
       << ", y1 = " << y1
       << ", x2 = " << x2
       << ", y2 = " << y2
       << ", x3 = " << x3
       << ", y3 = " << y3;
                     
     throw InvalidParameter( s.str() );
     
     
  }
  

   _R = 0;
   _centerX = 0.;
   _centerY = 0.;


   _x1 = x1;
   _y1 = y1;
   _x2 = x2;
   _y2 = y2;
   _x3 = x3;
   _y3 = y3;
               
               
   //check if x1 and x2 or x2 and x3 are equal. If they are, swap them around, so that those are not 0. (or else the slopes get infinite)
   // note that x1==x2==x3 is not possible as they would need to be on a line for that (and parallel, which we checked)
   
   if ( x1 == x2 ) {  // x1 and x2 have the same x --> we swap the points around (still it stays the same circle) so the 
                        // that they are now x1 and x3. because the line x1->x3 (i.e. its slope) isn't used in the calculations --> we don't care if it's zero.
      

      _x2 = x3;
      _y2 = y3;
      _x3 = x2;
      _y3 = y2;
      
   }
   else if ( x2 == x3 ) {  // x1 and x2 have the same x --> we swap the points around (still it stays the same circle) so the 
                        // that they are now x1 and x3. because the line x1->x3 (i.e. its slope) isn't used in the calculations --> we don't care if it's zero.
      

      _x2 = x1;
      _y2 = y1;
      _x1 = x2;
      _y1 = y2;
      
   }
   
   
   
   double ma = (_y2-_y1)/(_x2-_x1); //slope
   double mb = (_y3-_y2)/(_x3-_x2);
               
               
   _centerX = ( ma*mb*(_y1-_y3) + mb*(_x1+_x2) - ma*(_x2+_x3) )/( 2.*(mb-ma));
   _centerY = (-1./ma) * ( _centerX - (_x1+_x2)/2. ) + (_y1+_y2)/2;
            
   _R = sqrt (( _x1 - _centerX )*( _x1 - _centerX ) + ( _y1 - _centerY )*( _y1 - _centerY ));



  
  
  
  
}

