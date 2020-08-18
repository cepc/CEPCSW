#ifndef Crit4_NoZigZag_h
#define Crit4_NoZigZag_h


#include "Criteria/ICriterion.h"

namespace KiTrack{
   
   /** Criterion: forbids zig zagging: measure the angles in the xy plane, transpose them to the range from -pi to pi
    * and multiply: if there is a zigzag, the sign of the angle switches and the product of both angles becomes
    * negative.
    */
   class Crit4_NoZigZag : public ICriterion{
      
      
      
   public:
      
      /**
       * @param prodMin the minimum product of the two angles in degrees
       * 
       * @param prodMax the maxinum product of the two angles in degrees
       */
      Crit4_NoZigZag ( float prodMin , float prodMax );
      
      virtual bool areCompatible( Segment* parent , Segment* child );
      
      virtual ~Crit4_NoZigZag(){};
      
   private:
      
      float _prodMin{};
      float _prodMax{};
      
   };
   
}

#endif
