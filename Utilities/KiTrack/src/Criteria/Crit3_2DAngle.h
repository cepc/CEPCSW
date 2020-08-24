#ifndef Crit3_2DAngle_h
#define Crit3_2DAngle_h


#include "Criteria/ICriterion.h"

namespace KiTrack{

   /** Criterion: the angle between two 2-segments in the xy - plane
    */
   class Crit3_2DAngle : public ICriterion{



   public:
      
      /**
       * @param angleMax the maximum angle between 2 2-segments in grad
       */
      Crit3_2DAngle ( float angleMin, float angleMax );
      
      virtual bool areCompatible( Segment* parent , Segment* child );

      virtual ~Crit3_2DAngle(){};
    
   private:
      
      float _cosAngleMin{};
      float _cosAngleMax{};
      
   };

}

#endif

