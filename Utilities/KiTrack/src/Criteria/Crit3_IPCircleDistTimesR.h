#ifndef Crit3_IPCircleDistTimesR_h
#define Crit3_IPCircleDistTimesR_h

#include "Criteria/ICriterion.h"


namespace KiTrack{
   
   /** Criterion: the distance of the circle formed by the two segments from the IP multiplied by R
    */
   class Crit3_IPCircleDistTimesR : public ICriterion{
      
      
      
   public:
      
      Crit3_IPCircleDistTimesR ( float distToCircleMin , float distToCircleMax );
      
      virtual bool areCompatible( Segment* parent , Segment* child );
      
      virtual ~Crit3_IPCircleDistTimesR(){};
      
      
   private:
      
      float _distToCircleMax{};
      float _distToCircleMin{};
      
      
   };
   
}






#endif

