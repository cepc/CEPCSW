#ifndef Crit4_DistToExtrapolation_h
#define Crit4_DistToExtrapolation_h


#include "Criteria/ICriterion.h"

namespace KiTrack{
   
   /** Criterion: use the first 3 hits to extrapolate the location in xy for a given z of the last hit.
    * Then measure the distance from the extrapolation to the hit. (Also divide by the z distance to the last hit,
    * in order to take into account that with farther distances the accuracy drops)
    */
   class Crit4_DistToExtrapolation : public ICriterion{
      
      
      
   public:
      
      /**
       * @param distMax 
       */
      Crit4_DistToExtrapolation ( float distMin , float distMax );
      
      virtual bool areCompatible( Segment* parent , Segment* child );
      
      virtual ~Crit4_DistToExtrapolation(){};
      
   private:
      
      float _distMin{};
      float _distMax{};
      
   };
   
}

#endif
