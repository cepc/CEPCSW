#ifndef Crit2_DeltaTheta_MV_h
#define Crit2_DeltaTheta_MV_h


#include "Criteria/ICriterion.h"

namespace KiTrack{

   /** Criterion: the difference between the \f$ \theta \f$ angles of two hits in degrees. 
    * The \f$ \theta \f$ is the angle in the xy plane w.r.t. the positive x axis.
    * \f[ \theta = atan2(y,x) \f]
    */
   class Crit2_DeltaTheta_MV : public ICriterion{



   public:
      
      Crit2_DeltaTheta_MV ( float deltaThetaMin , float deltaThetaMax );
      
      virtual bool areCompatible( Segment* parent , Segment* child );

      virtual ~Crit2_DeltaTheta_MV(){};

    
   private:
      
      float _deltaThetaMax{};
      float _deltaThetaMin{};
      
      
      
   };

}

#endif

