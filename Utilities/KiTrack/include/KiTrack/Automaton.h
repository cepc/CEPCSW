#ifndef Automaton_h
#define Automaton_h

#include <vector>
#include "KiTrack/Segment.h"
#include "Criteria/ICriterion.h"



namespace KiTrack{
   
   
   /** A class for the cellular automaton.
    * 
    * For detailed info on the Cellular Automaton and its use for track reconstruction see 
    * <a href="../CellularAutomaton.pdf">Introduction to the Cellular Automaton</a>
    * 
    * <h3>About the class </h3>
    * 
    * This Cellular Automaton is specifically designed for track reconstruction, so it differs in the functionality
    * from other Cellular Automata, like those used for simulation of biological cells.
    * 
    * The basic information on how the CA works and what it needs is describet in the pdf above. Here I'll only sum up
    * the basics of this class.
    * 
    * The cells the Cellular Automaton deals with are here called segments as we deal with track reconstruction.
    * These are the main information entities the CA deals with and therefore the most important member variable is 
    * a container of these segments. (For more details on the segments see the doxygen of the Segment class)
    * To sum it up: Segments consist of hits. So they are a bit like tracks. Or if they only consist of one single hit
    * (so called 1-hit-segments) then they are like hits. 
    * But there a few features that distinguish them from hits or tracks
    *   - They can have parents and children. These are simply segments on top or below them. The idea is, that all
    * segments on the inside are stored as children if they could belong to the same track. And the same for all segments
    * on the outside possible belonging to the same track (the parents). So if we start at a segment and go to a child
    * and on to a grandchild and so on, we follow a possible track. If we erase those connections we decrease the number
    * of possible tracks.
    *   - Segments must have a layer marking where it is. (layer 0: inside, higher layers further outside). This is a
    * feature needed by the Automaton. Because it is an algorithm with discrete entities, we need some sort of discrete 
    * ordering. (If we for example had a Cellular Automaton for simulation of 2 dimensional cells and we arranged the 
    * cells for example on a chessboard, we always know the next cells. Here we have a 1 dimensional situation so we use
    * layers)
    *   - Segments have states: this is simply an integer number (an unsigned to be more precise). It is needed by the 
    * Automaton to find connections that go all the way through (see pdf!)
    *
    * The Segments can be added via the addSegment() method and are stored layerwise.
    * 
    * Once the Segments are all stored in the Cellular Automaton it can perform.
    * Via the method doAutomaton() it raises the states of the Segments until no change happens anymore.
    * When this is done Segments not connected all the way through can be discarded by the method 
    * cleanBadStates(). This reduces the number of possible tracks.
    * 
    * To get an initial Cellular Automaton to start with, the class SegmentBuilder can be used.
    * (It takes hits builds segments from them and establishes the first parent-child relations)
    * 
    * In order to sort out even more it is possible to go to longer segments. So instead of checking 1-hit-segments,
    * we can have a look at 2-hit-segments or 3-hit-segments. (And sort out much more along the way)
    * For this the method lengthenSegments() is used. It combines connected segments and creates segments from them,
    * that are exactly one hit longer. And then it is again time for connecting them.
    * (I distinguish here between connecting: "storing the link" and combining: "making 1 new segment out of two others")
    * 
    * When segments are connected only connections that make sense are made. This is really important! If we don't make
    * assumptions here, what Segments could belong together (i.e. could form a sensible track) we get lost in combinatorics.
    * For this the so called Criteria are used. Via the method addCriterion() (or addCriteria() ) they can be added
    * to the Cellular Automaton. All a Criterion does (see the Criterion doxygen for more info) is to say whether two
    * Segments would give a good match. These criteria can be anything that makes sense in distinguishing between 
    * real tracks and combinatorial background. For example if the segments are already long enough (already little tracks)
    * one can compare the radius of their helices or the angle under which they meet and so on. 
    * 
    * Once we have longer segments we can again use the doAutomaton() method and then cleanBadStates(). Then we lenghten
    * them once more and so on.
    * 
    * Between two such runs, one should of course clear the old Criteria with clearCriteria and reset the states of the
    * Segments with resetStates().
    *
    * That's it. In the end, when nothing more is to be done in the Cellular Automaton, we need to extract the track
    * candidates, it found. This is done via the getTracks() method. 
    *
    */ 
   class Automaton{
      
      
     
   public:
      
      Automaton(): _nConnections(0){}
      
      /**
       * delete all the segments
       */
      ~Automaton();
      
      
      /** Adds a segment to the automaton.\
       * Take care to set the layer of the segment before adding!
       */
      void addSegment ( Segment* segment );
      
      /**Lengthens the segments by one via adding the first hit of the next segment it is connected to
       * to it.
       * Also connects those longer segments with each other. ( one becomes a parent and one a child )
       * Segments that don't have connected segments to use to get longer, will die here. 
       */
      void lengthenSegments();
      
      /**Adds a criteria to the automaton. So it will be used, when the methods doAutomaton()
       * or cleanBadConnections() are called.
       */
      void addCriterion ( ICriterion* criterion ){ _criteria.push_back( criterion ); }
      void addCriteria ( std::vector< ICriterion* > criteria ){ _criteria.insert( _criteria.end() , criteria.begin() , criteria.end() ); }
      void clearCriteria() { _criteria.clear(); };
      
      /** Does iteration until the states of the segments don't change anymore.
       * 
       * In one iteration all segments are checked, if they have a neighbor.
       * A neighbor:
       * - has the same state
       * - fulfills all the criteria (this is automatically provided, because in lengthenSegments only connections that fullfill the Criteria are made)
       * 
       * If it has a neighbor its state will be raised by one at the end of the iteration.
       * (At the end only in theory: in the program it will be during the iteration, but
       * in a way, that it doesn't affect other segments)
       * 
       * When after an iteration the states didn't change, it stops.
       */
      void doAutomaton();
      
      /**
       * Erases all segments that don't have a state corresponding to their layer.
       * 
       * After the automaton is performed every segment, that has a neighbor, that has a neighbor, that...
       * ...that reaches layer 0 should have a state equal to its layer. All the others are not connected
       * to layer 0. All those get deleted by this method.
       * (Of course all the connections to them are erased as well)
       */
      void cleanBadStates();
      
      /**
       * Erase alls connections between segments, that don't satisfy the criteria.
       */
      void cleanBadConnections();
      
      /**Resets all the states of the segmens to 0 by calling the resetState() method of the segment
       * Also sets all segments back to active.
       */
      void resetStates();
      
      
      
      /** Get all the possible tracks in the automaton.
       * 
       * Tracks are built by starting from segments which don't have a parent (and are therefore the begin of a track).
       * For all children they have a new own track is created. The children themselves have children again, so the 
       * tracks split up again.
       * If we have for example a segment with 3 children, which each have 5 children, which each have 7 children,
       * we will get 1*3*5*7=105 tracks.
       * 
       * @return All tracks that are possible with the given segments and their connections. Tracks are returned
       * as a vector of hits. So the output will be a vector of a vector of hits
       * 
       * @param minHits the minimum number of hits that a track needs to have. All possible tracks,
       * that have less won't be considered as tracks and won't be returned.
       * 
       */
      //std::vector < std::vector< IHit* > > getTracks( unsigned minHits = 3 );
      std::vector < std::vector< IHit* > > getTracks( unsigned minHits = 2 ); // YV, 2 mini-vector hits can form a track     
      
      /**Returns all the tracks starting from this segment.
       * It is a recursive method and gets invoked by getTracks.
       */
      //std::vector < std::vector< IHit* > > getTracksOfSegment ( Segment* segment, std::vector< IHit* > hits , unsigned minHits = 3 );
      std::vector < std::vector< IHit* > > getTracksOfSegment ( Segment* segment, std::vector< IHit* > hits , unsigned minHits = 2 );    // YV, 2 mini-vector hits can form a track   
      
      /**
       * @return All the segments currently saved in the automaton
       */
      std::vector <const Segment*> getSegments() const;
      
      unsigned getNumberOfConnections(){ return _nConnections; }
      
   private:
      
      /** Here the segments are stored.
       * The vector corresponds to the layer.
       * The list corresponds to the segments on the layer.
       * _segments[2] is a list with all segments on layer 2.
       * 
       * The segments will be deleted by the Automaton in the destructor
       * 
       */
      std::vector < std::list < Segment* > > _segments{};
      
      /** A vector containing all the criteria, that are used in the Automaton
       */
      std::vector < ICriterion* > _criteria{};
      
      unsigned _nConnections{};
      
      
      
   };  
   
   
   
   
   
}






#endif


