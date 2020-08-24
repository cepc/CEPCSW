#ifndef Crit4_RChange_h
#define Crit4_RChange_h


#include "Criteria/ICriterion.h"

namespace KiTrack{
   
   /** Criterion: check for the change of the radii of the circles that can be made from the 3-hit segments
    */
   class Crit4_RChange : public ICriterion{
      
      
      
   public:
      
      /**
       * @param changeMax 
       */
      Crit4_RChange ( float changeMin , float changeMax );
      
      virtual bool areCompatible( Segment* parent , Segment* child );
      
      virtual ~Crit4_RChange(){};
      
   private:
      
      float _changeMin{};
      float _changeMax{};
      
   };
   
}

#endif
