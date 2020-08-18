#ifndef Crit3_3DAngle_h
#define Crit3_3DAngle_h


#include "Criteria/ICriterion.h"

namespace KiTrack{

   /** Criterion: the angle between two 2-segments
    */
   class Crit3_3DAngle : public ICriterion{



   public:
      
      /**
       * @param angleMax the maximum angle between 2 2-segments in grad
       */
      Crit3_3DAngle ( float angleMin, float angleMax );
      
      virtual bool areCompatible( Segment* parent , Segment* child );

      virtual ~Crit3_3DAngle(){};
    
   private:
      
      float _cosAngleMin{};
      float _cosAngleMax{};
      
   };

}

#endif

