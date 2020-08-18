#ifndef Crit2_RZRatio_h
#define Crit2_RZRatio_h


#include "Criteria/ICriterion.h"

namespace KiTrack{

   /** Criterion: distance of two hits divided by their z-distance.
    * \f[ \frac{\sqrt{ \Delta x^2 + \Delta y^2 + \Delta z^2 }}{\left| \Delta z \right|}\f]
    */
   class Crit2_RZRatio : public ICriterion{



   public:
      
      Crit2_RZRatio ( float ratioMin, float ratioMax );
      
      virtual bool areCompatible( Segment* parent , Segment* child );
      
      virtual ~Crit2_RZRatio(){};
    
   private:
      
      float _ratioMax{};
      float _ratioMin{};
      
      
   };

}

#endif

