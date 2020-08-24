#ifndef Crit2_PT_MV_h
#define Crit2_PT_MV_h

#include "Criteria/ICriterion.h"


namespace KiTrack{
   
   /** Criterion: the transversal momentum 
    */
   class Crit3_PT_MV : public ICriterion{
      
      
      
   public:
      
      Crit3_PT_MV ( float ptMin , float ptMax , float Bz = 3.5 );
      
      virtual bool areCompatible( Segment* parent , Segment* child );
      
      virtual ~Crit3_PT_MV(){};
      
      
   private:
      
      float _ptMin{};
      float _ptMax{};
      float _Bz{};
      
      
   };
   
}













#endif

