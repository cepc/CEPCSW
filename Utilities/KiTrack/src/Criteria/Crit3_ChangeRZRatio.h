#ifndef Crit2_ChangeRZRatio_h
#define Crit2_ChangeRZRatio_h

#include "Criteria/ICriterion.h"


namespace KiTrack{
   
   /** Criterion: the change of the ratio of the distance of two hits over the z distance
    */
   class Crit3_ChangeRZRatio : public ICriterion{
      
      
      
   public:
      
      Crit3_ChangeRZRatio ( float minChange , float maxChange );
      
      virtual bool areCompatible( Segment* parent , Segment* child );
      
      virtual ~Crit3_ChangeRZRatio(){};
      
      
   private:
      
      float _ratioChangeMaxSquared{};
      float _ratioChangeMinSquared{};

      
      
   };
   
}













#endif

