#ifndef Crit2_HelixWithIP_h
#define Crit2_HelixWithIP_h


#include "Criteria/ICriterion.h"

namespace KiTrack{

   /** Criterion: Check if two hits are compatible with a helix crossing the IP.
    * The procedure is this: a circle is calculated from the two hits and the IP.
    * The angle between the center of the circle and two hits is calles alpha here.
    * alpha between IP and the first hit is proportional to the z-distance of the IP
    * and the first hit in the same way alpha between the first and second hit should 
    * be proportional to the z-distance of those hits. The value calculated is
    * \f[ \frac{ \frac{\alpha_1}{\Delta z_1} }{ \frac{\alpha_2}{\Delta z_2} } \simeq 1 \f]
    * and is 1 for a perfect helix around the z - axis.
    */
   class Crit2_HelixWithIP : public ICriterion{



   public:
      
      Crit2_HelixWithIP ( float ratioMin , float ratioMax );
      
      virtual bool areCompatible( Segment* parent , Segment* child );

      virtual ~Crit2_HelixWithIP(){};

    
   private:
      
      float _ratioMax{};
      float _ratioMin{};
      
      
      
   };

}

#endif

