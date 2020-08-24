#ifndef Crit2_PT_h
#define Crit2_PT_h

#include "Criteria/ICriterion.h"


namespace KiTrack{
   
   /** Criterion: the transversal momentum 
    */
   class Crit3_PT : public ICriterion{
      
      
      
   public:
      
      Crit3_PT ( float ptMin , float ptMax , float Bz = 3.5 );
      
      virtual bool areCompatible( Segment* parent , Segment* child );
      
      virtual ~Crit3_PT(){};
      
      
   private:
      
      float _ptMin{};
      float _ptMax{};
      float _Bz{};
      
      
   };
   
}













#endif

