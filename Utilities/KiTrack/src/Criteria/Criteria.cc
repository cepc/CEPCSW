#include "Criteria/Criteria.h"
#include "Criteria/ICriterion.h"

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

using namespace KiTrack;

std::set< std::string > Criteria::getAllCriteriaNames(){

   
   std::set< std::string > critNames;
   
   critNames.insert( "Crit2_RZRatio" );
   critNames.insert( "Crit2_StraightTrackRatio" );
   critNames.insert( "Crit2_DeltaPhi" );
   critNames.insert( "Crit2_HelixWithIP" );
   critNames.insert( "Crit2_DeltaRho" );

   critNames.insert( "Crit3_ChangeRZRatio" );
   critNames.insert( "Crit3_PT" );
   critNames.insert( "Crit3_2DAngle" );
   critNames.insert( "Crit3_2DAngleTimesR" );
   critNames.insert( "Crit3_3DAngle" );
   critNames.insert( "Crit3_3DAngleTimesR" );
   critNames.insert( "Crit3_IPCircleDist" );
   critNames.insert( "Crit3_IPCircleDistTimesR" );
   

   critNames.insert( "Crit4_2DAngleChange" );
   critNames.insert( "Crit4_3DAngleChange" );
   critNames.insert( "Crit4_3DAngleChangeNormed" );
   critNames.insert( "Crit4_DistToExtrapolation" );
   critNames.insert( "Crit4_PhiZRatioChange" );
   critNames.insert( "Crit4_DistOfCircleCenters" );
   critNames.insert( "Crit4_NoZigZag" );
   critNames.insert( "Crit4_RChange" );

   // MiniVector based Cellular Automaton for VXD

   critNames.insert( "Crit2_DeltaPhi_MV" );
   critNames.insert( "Crit2_Distance_MV" );
   critNames.insert( "Crit2_DeltaTheta_MV" ); 
   critNames.insert( "Crit3_NoZigZag_MV" );
   critNames.insert( "Crit3_PT_MV" );

   return critNames;

}


std::set < std::string > Criteria::getTypes(){
 
   std::set< std::string > critNames = getAllCriteriaNames();
   std::set< std::string > types;
   
   std::set< std::string >::iterator it;
   
   
   for( it = critNames.begin(); it != critNames.end(); it++ ){
      
      
      ICriterion* crit = Criteria::createCriterion( *it );
      
      types.insert( crit->getType() );
      
      delete crit;
      
      
   }
   
   return types;
   
   
}


std::set< std::string > Criteria::getCriteriaNames( std::string type ){
   
   
   std::set< std::string > criteria;
   std::set< std::string > critNames = getAllCriteriaNames();
   
   
   std::set< std::string >::iterator it;
   
   for( it = critNames.begin(); it != critNames.end(); it++ ){
      
      
      ICriterion* crit = Criteria::createCriterion( *it );
      
      if ( crit->getType() == type ) criteria.insert( *it );
      
      delete crit;
      
   }

   return criteria;
      
    
   
   
}


ICriterion* Criteria::createCriterion( std::string critName, float min , float max ) {
   
   
   
   if ( critName == "Crit2_RZRatio" ) return ( new Crit2_RZRatio( min , max ) );
   
   else if ( critName == "Crit2_StraightTrackRatio" ) return ( new Crit2_StraightTrackRatio( min , max ) );
   
   else if ( critName == "Crit2_DeltaPhi" ) return ( new Crit2_DeltaPhi( min , max ) );
   
   else if ( critName == "Crit2_HelixWithIP" ) return ( new Crit2_HelixWithIP( min , max ) );
   
   else if ( critName == "Crit2_DeltaRho" ) return ( new Crit2_DeltaRho( min , max ) );
   
   else if ( critName == "Crit3_ChangeRZRatio" ) return ( new Crit3_ChangeRZRatio( min , max ) );
   
   else if ( critName == "Crit3_PT" ) return ( new Crit3_PT( min , max ) );
   
   else if ( critName == "Crit3_2DAngle" ) return ( new Crit3_2DAngle( min , max ) );
   
   else if ( critName == "Crit3_2DAngleTimesR" ) return ( new Crit3_2DAngleTimesR( min , max ) );
   
   else if ( critName == "Crit3_3DAngle" ) return ( new Crit3_3DAngle( min , max ) );
   
   else if ( critName == "Crit3_3DAngleTimesR" ) return ( new Crit3_3DAngleTimesR( min , max ) );
   
   else if ( critName == "Crit3_IPCircleDist" ) return ( new Crit3_IPCircleDist( min , max ) );
   
   else if ( critName == "Crit3_IPCircleDistTimesR" ) return ( new Crit3_IPCircleDistTimesR( min , max ) );
   
   else if ( critName == "Crit4_2DAngleChange" ) return ( new Crit4_2DAngleChange( min , max ) );
   
   else if ( critName == "Crit4_3DAngleChange" ) return ( new Crit4_3DAngleChange( min , max ) );
   
   else if ( critName == "Crit4_3DAngleChangeNormed" ) return ( new Crit4_3DAngleChangeNormed( min , max ) );
   
   else if ( critName == "Crit4_DistToExtrapolation" ) return ( new Crit4_DistToExtrapolation( min , max ) );
   
   else if ( critName == "Crit4_PhiZRatioChange" ) return ( new Crit4_PhiZRatioChange( min , max ) );
   
   else if ( critName == "Crit4_DistOfCircleCenters" ) return ( new Crit4_DistOfCircleCenters( min , max ) );
   
   else if ( critName == "Crit4_NoZigZag" ) return ( new Crit4_NoZigZag( min , max ) );
   
   else if ( critName == "Crit4_RChange" ) return ( new Crit4_RChange( min , max ) );

   // Mini-Vector based

   else if ( critName == "Crit2_DeltaPhi_MV" ) return ( new Crit2_DeltaPhi_MV( min , max ) );

   else if ( critName == "Crit2_Distance_MV" ) return ( new Crit2_Distance_MV( min , max ) );

   else if ( critName == "Crit2_DeltaTheta_MV" ) return ( new Crit2_DeltaTheta_MV( min , max ) );

   else if ( critName == "Crit3_NoZigZag_MV" ) return ( new Crit3_NoZigZag_MV( min , max ) );

   else if ( critName == "Crit3_PT_MV" ) return ( new Crit3_PT_MV( min , max ) );
   
   else {
      
      std::string s = "Criteria::The criterion \"" + critName + 
                      "\" is not known. Make sure the class Criteria has this criterion listed in the createCriterion method";
      
      throw UnknownCriterion( s );
      
   }

      
      
    
}



std::vector< std::string > Criteria::getAllCriteriaNamesVec(){
   
   std::vector < std::string > allCriteriaNamesVec;
   std::set< std::string > critNames = getAllCriteriaNames();
   
   
   std::set< std::string >::iterator it;
   
   for( it = critNames.begin(); it != critNames.end(); it++ ){
   
      
      allCriteriaNamesVec.push_back( *it );
      
   }
   
   return allCriteriaNamesVec;
   
}


void Criteria::getLeftRight( std::string critName, float & left, float & right ){
   
   
   if ( critName == "Crit2_RZRatio" ) { left = 0.5; right = 0.5; }
   
   else if ( critName == "Crit2_RZRatio" ) { left = 0.5; right = 0.5; }
   
   else if ( critName == "Crit2_StraightTrackRatio" ) { left = 0.5; right = 0.5; }
   
   else if ( critName == "Crit2_DeltaPhi" ) { left = 0.; right = 1.; }
   
   else if ( critName == "Crit2_HelixWithIP" ) { left = 0.5; right = 0.5; }
   
   else if ( critName == "Crit2_DeltaRho" ) { left = 0.5; right = 0.5; }
   
   else if ( critName == "Crit3_ChangeRZRatio" ) { left = 0.5; right = 0.5; }
   
   else if ( critName == "Crit3_PT" ) { left = 0.5; right = 0.5; }
   
   else if ( critName == "Crit3_2DAngle" ) { left = 0.; right = 1.; }
   
   else if ( critName == "Crit3_2DAngleTimesR" ) { left = 0.; right = 1.; }
   
   else if ( critName == "Crit3_3DAngle" ) { left = 0.; right = 1.; }
   
   else if ( critName == "Crit3_3DAngleTimesR" ) { left = 0.; right = 1.; }
   
   else if ( critName == "Crit3_IPCircleDist" ) { left = 0.; right = 1.; }
   
   else if ( critName == "Crit3_IPCircleDistTimesR" ) { left = 0.; right = 1.; }
   
   else if ( critName == "Crit4_2DAngleChange" ) { left = 0.5; right = 0.5; }
   
   else if ( critName == "Crit4_3DAngleChange" ) { left = 0.5; right = 0.5; }
   
   else if ( critName == "Crit4_3DAngleChangeNormed" ) { left = 0.5; right = 0.5; }
   
   else if ( critName == "Crit4_DistToExtrapolation" ) { left = 0.; right = 1.; }
   
   else if ( critName == "Crit4_PhiZRatioChange" ) { left = 0.5; right = 0.5; }
   
   else if ( critName == "Crit4_DistOfCircleCenters" ) { left = 0.; right = 1.; }
   
   else if ( critName == "Crit4_NoZigZag" ) { left = 0.5; right = 0.5; }
   
   else if ( critName == "Crit4_RChange" ) { left = 0.5; right = 0.5; }
   
   else { left = 0.5; right = 0.5; }
   
 
}
   

