#ifndef Crit4_2DAngleChange_h
#define Crit4_2DAngleChange_h


#include "Criteria/ICriterion.h"

namespace KiTrack{
   
   /** Criterion: change of the angle (in the xy plane) between segments in the xy plane
    */
   class Crit4_2DAngleChange : public ICriterion{
      
      
      
   public:
      
      /**
       * @param changeMax 
       */
      Crit4_2DAngleChange ( float changeMin , float changeMax );
      
      virtual bool areCompatible( Segment* parent , Segment* child );
      
      virtual ~Crit4_2DAngleChange(){};
      
   private:
      
      float _changeMin{};
      float _changeMax{};
      
   };
   
}

#endif
