#include "Tools/KiTrackMarlinCEDTools.h"

#include <cmath>

#include <CLHEP/Random/RandFlat.h>

#include "MarlinCED.h"


void KiTrackMarlin::drawAutomatonSegments( const Automaton& automaton ){
   
   
   std::vector< const Segment* > segments = automaton.getSegments();
   
   
   for( unsigned i=0 ; i < segments.size(); i++ ){
      
      
      const Segment* segment = segments[i];
      std::vector < IHit* > hits = segment->getHits();         
      
     
      if ( hits.size() == 1){ //exactly one hit, so draw a point
         
         
         IHit* a = hits[0];
         ced_hit( a->getX() ,a->getY() , a->getZ() , 0 , 3 ,0xff0000 );
         
         
      }
      else{ //more than one point or no points
         
         for( unsigned j=1 ; j< hits.size() ; j++ ){ // over all hits in the segment (as we connect it with the previous we start with hit 1)
            
            IHit* a = hits[j];
            IHit* b = hits[j-1];
            
            
            unsigned int color=0;
            unsigned int red=0;
            unsigned int blue=0;
            unsigned int green=0;
            
            float p =  sqrt ((float)  segment->getInnerState() / (float) ( 7 ));
            
            green = unsigned( ceil ( (1.-p) * 255 ) );
            red = unsigned( floor( 255*p ) );
            blue = unsigned( ceil ( (1.-p) * 255 ) );
            
            color = red * 256*256 + green * 256 + blue;
            
            
            ced_line_ID( a->getX() ,a->getY() , a->getZ() , b->getX() ,b->getY() , b->getZ() , 2 , segment->getInnerState()+2 , color, 0);
            
         }
         
      }
      
   }
   
   
}



void KiTrackMarlin::drawTrack( ITrack* track, int color ){
   
   
   std::vector < IHit* > hits = track->getHits();     
   
   for( unsigned j=1 ; j< hits.size() ; j++ ){ // over all hits in the segment (as we connect it with the previous we start with hit 1)
      
      IHit* a = hits[j];
      IHit* b = hits[j-1];
      
      ced_line_ID( a->getX() ,a->getY() , a->getZ() , b->getX() ,b->getY() , b->getZ() , 2 , 5 , color, 0);
      
   }
   
   
   
   
}


void KiTrackMarlin::drawTrackRandColor( ITrack* track ){
   
   
   int red = 0;
   int green = 0;
   int blue = 0;
   
   int dist = 0;
   int distmin = 60;
   
   while( dist < distmin ){
      
      red = CLHEP::RandFlat::shootInt(256);
      green = CLHEP::RandFlat::shootInt(256);
      blue = CLHEP::RandFlat::shootInt(256);
      
      
      dist = std::abs( red - green );
      if( std::abs( green - blue ) > dist ) dist = std::abs( green - blue );
      if( std::abs( red - blue ) > dist ) dist = std::abs( red - blue );
      
   }
   
   int color = red * 256*256 + green * 256 + blue;
   drawTrack( track, color );
   
   
   
}








