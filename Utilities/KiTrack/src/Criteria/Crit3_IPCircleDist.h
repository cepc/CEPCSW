#ifndef Crit2_IPCircleDist_h
#define Crit2_IPCircleDist_h

#include "Criteria/ICriterion.h"


namespace KiTrack{
   
   /** Criterion: the distance of the IP from the circle the 3 hits form (in the xy plane)
    */
   class Crit3_IPCircleDist : public ICriterion{
      
      
      
   public:
      
      Crit3_IPCircleDist ( float distToCircleMin , float distToCircleMax );
      
      virtual bool areCompatible( Segment* parent , Segment* child );
      
      virtual ~Crit3_IPCircleDist(){};
      
      
   private:
      
      float _distToCircleMax{};
      float _distToCircleMin{};
      
      
   };
   
}













#endif

