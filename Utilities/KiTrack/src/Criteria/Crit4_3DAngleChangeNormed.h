#ifndef Crit4_3DAngleChangeNormed_h
#define Crit4_3DAngleChangeNormed_h


#include "Criteria/ICriterion.h"

namespace KiTrack{
   
   /** Criterion: check for the change of the 3D angle and normalise it with R
    */
   class Crit4_3DAngleChangeNormed : public ICriterion{
      
      
      
   public:
      
      /**
       * @param changeMax 
       */
      Crit4_3DAngleChangeNormed ( float changeMin , float changeMax );
      
      virtual bool areCompatible( Segment* parent , Segment* child );
      
      virtual ~Crit4_3DAngleChangeNormed(){};
      
   private:
      
      float _changeMin{};
      float _changeMax{};
      
   };
   
}

#endif
