#ifndef KiTrackMarlinCEDTools_h
#define KiTrackMarlinCEDTools_h

#include "KiTrack/Automaton.h"
#include "KiTrack/ITrack.h"

using namespace KiTrack;

namespace KiTrackMarlin{
   
   
   
   /**
   * Draws the segments of an automaton.
   * Segments will be colored according to their state.
   * Also the higher the state, the thicker the line.
   */
   void drawAutomatonSegments( const Automaton& automaton );
   
   /**
    * Draws a track by making straight lines between the hits
    */
   void drawTrack( ITrack* track, int color = 0x00ff00 );
   
   /** Draw a track in a random color, but not too bright */
   void drawTrackRandColor( ITrack* track );
   
}

#endif 

