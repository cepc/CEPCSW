#ifndef Crit3_3DAngleTimesR_h
#define Crit3_3DAngleTimesR_h


#include "Criteria/ICriterion.h"

namespace KiTrack{

   /** Criterion: the angle between two 2-segments multiplied by the radius of the circle the segments form
    */
   class Crit3_3DAngleTimesR : public ICriterion{



   public:
      
      /**
       * @param angleMax the maximum angle between 2 2-segments in grad times the radius of the circle the segments form
       */
      Crit3_3DAngleTimesR ( float angleMin, float angleMax );
      
      virtual bool areCompatible( Segment* parent , Segment* child );

      virtual ~Crit3_3DAngleTimesR(){};
    
   private:
      
      float _angleMin{};
      float _angleMax{};
      
   };

}

#endif

