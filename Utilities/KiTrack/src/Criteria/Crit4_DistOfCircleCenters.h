#ifndef Crit4_DistOfCircleCenters_h
#define Crit4_DistOfCircleCenters_h


#include "Criteria/ICriterion.h"

namespace KiTrack{
   
   /** Criterion: make circles from the semgments and measure the distances of their centers
    */
   class Crit4_DistOfCircleCenters : public ICriterion{
      
      
      
   public:
      
      /**
       * @param distMax 
       */
      Crit4_DistOfCircleCenters ( float distMin , float distMax );
      
      virtual bool areCompatible( Segment* parent , Segment* child );
      
      virtual ~Crit4_DistOfCircleCenters(){};
      
   private:
      
      float _distMax{};
      float _distMin{};
      
   };
   
}

#endif
