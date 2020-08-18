#ifndef Crit3_2DAngleTimesR_h
#define Crit3_2DAngleTimesR_h


#include "Criteria/ICriterion.h"

namespace KiTrack{

   /** Criterion: the angle in the xy plane between two 2-segments multiplied by the radius of the circle they form
    */
   class Crit3_2DAngleTimesR : public ICriterion{



   public:
      
      /**
       * @param angleMax the maximum angle between 2 2-segments in grad times radius
       */
      Crit3_2DAngleTimesR ( float angleMin, float angleMax );
      
      virtual bool areCompatible( Segment* parent , Segment* child );

      virtual ~Crit3_2DAngleTimesR(){};
    
   private:
      
      float _angleMin{};
      float _angleMax{};
      
   };

}

#endif

