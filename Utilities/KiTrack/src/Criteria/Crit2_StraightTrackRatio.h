#ifndef Crit2_StraightTrackRatio_h
#define Crit2_StraightTrackRatio_h


#include "Criteria/ICriterion.h"

namespace KiTrack{

   /** Criterion: for straight tracks: if the line between the two hits points towards IP.
    * Calculated is
    * 
    * \f[ \frac{ \frac{\rho_1}{z_1} }{ \frac{\rho_2}{z_2} } \simeq 1 \f]
    * 
    * , where \f$ \rho = \sqrt{ x^2 + y^2 }\f$
    */
   class Crit2_StraightTrackRatio : public ICriterion{



   public:
      
      Crit2_StraightTrackRatio ( float ratioMin, float ratioMax );
      
      virtual bool areCompatible( Segment* parent , Segment* child );

      virtual ~Crit2_StraightTrackRatio(){};

    
   private:
      
      float _ratioMax{};
      float _ratioMin{};
      
      
      
   };

}

#endif

