#ifndef SegmentBuilder_h
#define SegmentBuilder_h

#include "Criteria/ICriterion.h"
#include "KiTrack/ISectorConnector.h"
#include "KiTrack/Automaton.h"

namespace KiTrack{


   /** This classe builds the Cellular Automaton from the hits
    * 
    * It can be used to take all the autHits stored in a map< int , vector < IHits > > and makes
    * 1-hit-segments out of them ( see the class Segment for more info on them ).
    * 
    * The created 1-segments then are connected.
    * 
    * For the rules of connecting criteria and hitConnectors can be added to the object:
    * 
    * - a hitConnector takes the sector of the segment ( for example a cellID0 or a layer number or an own code ) and returns
    * all the sectors we might connect to. So there we get the information like: "this segment can be connected
    * to layer 3 and 4, module 7,8,9 in forward direction".
    * 
    * - the criteria take two segments and return whether they are compatible or not. 
    * A criterion could check for anything, that is stored in the segments. 
    * For example: if the line formed from two 1-hit segments passes close by the IP might
    * be a criterion for very stiff tracks.
    * 
    * So the hitConnectors tell us were to look and the criteria whether to connect. If two 1-hit segments are found,
    * that are compatible, they will be connected. Connected means: The inner 1-segment will save the outer one as a parent
    * and the outer one will save the inner one as a child.
    * 
    * All this (except adding hitConnectors and Criteria) is done in the method get1SegAutomaton.
    * 
    * This method finally then returns an automaton (segments sorted by their layers), ready to be used.
    * 
    */   
   class SegmentBuilder{
      
     
   public: 
      
      /**
       * @param ftdRep the FTDRepresentation to take the autHits from
       */
      SegmentBuilder(  std::map< int , std::vector< IHit* > > map_sector_hits );
      
      /** Adds a criterion. 
       */
      void addCriterion ( ICriterion* criterion ){ _criteria.push_back( criterion );};
      
      /** Adds criteria
       */
      void addCriteria ( std::vector< ICriterion* > criteria){ _criteria.insert( _criteria.end(), criteria.begin() , criteria.end() ); }
      
      /** Adds a hitConnector
       */
      void addSectorConnector ( ISectorConnector* connector ){ _sectorConnectors.push_back( connector ); };
      
      /**
       * @return An automaton containing all the hits from the FTDRepresentation sorted now by layers. 
       * (attention: those are not necessarily the same layers as the FTD layers).
       * Also they are connected.
       */
      Automaton get1SegAutomaton();
      
      
   private:
      
      
      std::vector <ICriterion* > _criteria{};
      std::vector <ISectorConnector* > _sectorConnectors{};
      
      std::map< int , std::vector< IHit* > > _map_sector_hits{};
      
      
      
      
      
   };





}


#endif





