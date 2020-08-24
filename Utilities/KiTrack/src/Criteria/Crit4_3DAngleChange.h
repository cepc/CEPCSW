#ifndef Crit4_3DAngleChange_h
#define Crit4_3DAngleChange_h


#include "Criteria/ICriterion.h"

namespace KiTrack{
   
   /** Criterion: change of the angle between segments
    */
   class Crit4_3DAngleChange : public ICriterion{
      
      
      
   public:
      
      /**
       * @param changeMax 
       */
      Crit4_3DAngleChange ( float changeMin , float changeMax );
      
      virtual bool areCompatible( Segment* parent , Segment* child );
      
      virtual ~Crit4_3DAngleChange(){};
      
   private:
      
      float _changeMin{};
      float _changeMax{};
      
   };
   
}

#endif
