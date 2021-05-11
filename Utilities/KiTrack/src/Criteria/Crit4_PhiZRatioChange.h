#ifndef Crit4_PhiZRatioChange_h
#define Crit4_PhiZRatioChange_h


#include "Criteria/ICriterion.h"

namespace KiTrack{
   
   /** Criterion: check for the change of \f$ \frac{\phi}{Z} \f$. (\f$ \phi \f$ = angle in the circle formed by the hits)
    */
   class Crit4_PhiZRatioChange : public ICriterion{
      
      
      
   public:
      
      /**
       * @param changeMax 
       */
      Crit4_PhiZRatioChange ( float changeMin , float changeMax );
      
      virtual bool areCompatible( Segment* parent , Segment* child );
      
      virtual ~Crit4_PhiZRatioChange(){};
      
   private:
      
      float _changeMin{};
      float _changeMax{};
      
   };
   
}

#endif
