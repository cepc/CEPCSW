#ifndef Crit3_NoZigZag_MV_h
#define Crit3_NoZigZag_MV_h


#include "Criteria/ICriterion.h"

namespace KiTrack{
   
   /** Criterion: forbids zig zagging: measure the angles in the xy plane, transpose them to the range from -pi to pi
    * and multiply: if there is a zigzag, the sign of the angle switches and the product of both angles becomes
    * negative.
    */
   class Crit3_NoZigZag_MV : public ICriterion{
      
      
      
   public:
      
      /**
       * @param prodMin the minimum product of the two angles in degrees
       * 
       * @param prodMax the maxinum product of the two angles in degrees
       */
      Crit3_NoZigZag_MV ( float prodMin , float prodMax );
      
      virtual bool areCompatible( Segment* parent , Segment* child );
      
      virtual ~Crit3_NoZigZag_MV(){};
      
   private:
      
      float _prodMin{};
      float _prodMax{};
      
   };
   
}

#endif
