#ifndef Crit2_DeltaPhi_h
#define Crit2_DeltaPhi_h


#include "Criteria/ICriterion.h"

namespace KiTrack{

   /** Criterion: the difference between the \f$ \phi \f$ angles of two hits in degrees. 
    * The \f$ \phi \f$ is the angle in the xy plane w.r.t. the positive x axis.
    * \f[ \phi = atan2(y,x) \f]
    */
   class Crit2_DeltaPhi : public ICriterion{



   public:
      
      Crit2_DeltaPhi ( float deltaPhiMin , float deltaPhiMax );
      
      virtual bool areCompatible( Segment* parent , Segment* child );

      virtual ~Crit2_DeltaPhi(){};

    
   private:
      
      float _deltaPhiMax{};
      float _deltaPhiMin{};
      
      
      
   };

}

#endif

