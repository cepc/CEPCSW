#ifndef Crit2_Distance_MV_h
#define Crit2_Distance_MV_h


#include "Criteria/ICriterion.h"

namespace KiTrack{

   /** Criterion: the difference between the \f$ \phi \f$ angles of two hits in degrees. 
    * The \f$ \phi \f$ is the angle in the xy plane w.r.t. the positive x axis.
    * \f[ \phi = atan2(y,x) \f]
    */
   class Crit2_Distance_MV : public ICriterion{



   public:
      
      Crit2_Distance_MV ( float deltaPos2Min , float deltaPos2Max );
      
      virtual bool areCompatible( Segment* parent , Segment* child );

      virtual ~Crit2_Distance_MV(){};

    
   private:
      
      float _deltaPos2Max{};
      float _deltaPos2Min{};
      
      
      
   };

}

#endif

