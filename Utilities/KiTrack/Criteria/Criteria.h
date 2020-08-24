#ifndef Criteria_h
#define Criteria_h
/*
#include "Criteria/Crit2_RZRatio.h"
#include "Criteria/Crit2_StraightTrackRatio.h"
#include "Criteria/Crit2_DeltaPhi.h"
#include "Criteria/Crit2_HelixWithIP.h"
#include "Criteria/Crit2_DeltaRho.h"

#include "Criteria/Crit3_ChangeRZRatio.h"  
#include "Criteria/Crit3_PT.h"
#include "Criteria/Crit3_2DAngle.h"
#include "Criteria/Crit3_2DAngleTimesR.h"
#include "Criteria/Crit3_3DAngle.h"
#include "Criteria/Crit3_3DAngleTimesR.h"
#include "Criteria/Crit3_IPCircleDist.h"  
#include "Criteria/Crit3_IPCircleDistTimesR.h"

#include "Criteria/Crit4_2DAngleChange.h"    
#include "Criteria/Crit4_3DAngleChange.h" 
#include "Criteria/Crit4_3DAngleChangeNormed.h"
#include "Criteria/Crit4_DistToExtrapolation.h"  
#include "Criteria/Crit4_PhiZRatioChange.h"
#include "Criteria/Crit4_DistOfCircleCenters.h"
#include "Criteria/Crit4_NoZigZag.h"
#include "Criteria/Crit4_RChange.h"

// Criteria for Mini - Vector based Cellular Automaton for VXD
#include "Criteria/Crit2_DeltaPhi_MV.h"
#include "Criteria/Crit2_Distance_MV.h"
#include "Criteria/Crit2_DeltaTheta_MV.h"
#include "Criteria/Crit3_NoZigZag_MV.h"
#include "Criteria/Crit3_PT_MV.h"
*/

#include <vector>
#include <set>
#include <string>

namespace KiTrack{
  class ICriterion;
   /*
    * Information about all Criteria.
    * 
    * For example bundles the includes.
    * 
    * Author: Robin Glattauer, HEPHY
    */   
   class Criteria {
      
      
      
   public:
      
      /** @return a vector of strings that represent all types of criteria stored.
       * For example: "2Hit" or "3Hit" or "4Hit"
       */
      static std::set< std::string > getTypes();
      
            
      /** @return a vector of all Criteria of a certain type
       * 
       * @param type the type of Criteria, that is wanted (e.g. "2Hit")
       */
      static std::set< std::string > getCriteriaNames( std::string type );
      
      /** @return the names of all Criteria in a set
       */
      static std::set< std::string > getAllCriteriaNames();
      
      /** A convenience method to get all the criteria in a vector (gives the same result as getAllCriteriaNames, but
       * instead of a set, returns it as a vector)
       * 
       * @return the names of all Criteria in a vector
       */
      static std::vector< std::string > getAllCriteriaNamesVec();
      
      
      /**
       * Creates a Criterion with the name and the min and max values
       * 
       * @return a "new" Criterion (i.e. needs to be deleted later on)
       */
      static ICriterion* createCriterion( std::string critName , float min=0. , float max=0. ) ;
      
      /**
       * Sets values for the passed referneced floats left and right. They indicate how
       * the specified criterion should be cut, if necessary. Say you want for example
       * a 99% quantile, so that 99% of your true tracks are within it.
       * A criterion like the angle between two segments then needs to define a boarder like:
       * between an angle of 1° and of 9° there will be 99%.
       * So 1% is outside. But should 1% be the ones with a bigger angle or with a smaller angle or should
       * this be 50:50? 
       * 
       * This is defined by left and right. Left is the proportion, that is taken away on the left side and
       * right is the one that is taken away on the right side.
       * In the case of an angle we will most probably have lots around 0° and a long tail to the right, so
       * left = 0 and right = 1 seems like a good idea.
       * Standard is of course 0.5 and 0.5
       */
      static void getLeftRight( std::string critName, float & left, float & right );
      
   
      
   };

}

#endif

