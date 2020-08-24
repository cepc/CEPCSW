#ifndef SimpleCircle_h
#define SimpleCircle_h

#include "KiTrack/KiTrackExceptions.h"

namespace KiTrack{

/** A simple class representing a circle.
 * 
 * In the constructor it builds a cricle from 3 2-dimensional points.
 * After that, the parameters can be read out (radius and position of the center)
 */
class SimpleCircle {
  
  
  
  public:
      
    
      SimpleCircle ( double x1 , double y1 , double x2 , double y2 , double x3, double y3 ) ;
      
      double getRadius() {return _R;};
      double getCenterX() {return _centerX;};
      double getCenterY() {return _centerY;};
      
  private:
      
      
      double _R{};
      double _centerX{};
      double _centerY{};
      
      double _x1{};
      double _x2{};
      double _x3{};
      
      double _y1{};
      double _y2{};
      double _y3{};



};




} // end of namespace KiTrack




#endif


