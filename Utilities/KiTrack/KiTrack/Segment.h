#ifndef Segment_h
#define Segment_h

#include <vector>
#include <list>
#include <string>

#include "KiTrack/IHit.h"

namespace KiTrack{

   
   /** A Segment is something like a track or a part of a track: it consists of hits linked together.
    * 
    * The simplest Segment consists of only one single hit (the truly smallest part of a track),
    * a 1-hit-segment.
    * A segment of two hits (a 2-hit-segment) is the next bigger version and so on.
    * 
    * Segments are used by the Cellular Automaton (see class Automaton), which uses them to find tracks.
    * 
    * The main difference to a hit (in case of 1-hit-segments) or a track (in case of segments with more hits) is, that
    * the segments can have connection to other Segments. They can have children and parents.
    * Children are connected Segments on the inside, Parents are connected Segments on the outside.
    * 
    * Inside and outside are w.r.t. the layer a segment is on. Every Segment has a layer (getLayer(), setLayer() ). The 
    * layer indicates the place of the segment (whereever that place is. e.g. a detector ). Layer 0 usually means inside 
    * and higher layers are further outside.
    * 
    * Also every Segment has a state. It is another feature requested by the Cellular Automaton. It gives the CA the possibility
    * to check the quality of a Segment and manipulate it (the higher the state, usually the better; but for details on
    * the states in the Cellular Autoamton see the <a href="../CellularAutomaton.pdf">Introduction to the Cellular Automaton</a> )
    * 
    * Segments can skip layers (a 1-hit-segment can't, since it is just a hit and can only be on one layer, but a 
    * 2-hit-segment could already be the connection between layer 2 and layer 4, thus skipping a layer).
    * For this the state is not only an integer, but a vector of integers representing states for every layer.
    * (the 2-hit-segment skipping layer number 3 can therefore be imagined as two 2-hit-segments, each of them having
    * a seperate state )
    * 
    */   
   class Segment {
   
   
   public:
      
      Segment( std::vector <IHit*> hits);
      Segment( IHit* hit);
      
      
      void deleteParent ( Segment* delParent ){ _parents.remove( delParent );};
      void deleteChild ( Segment* delChild ){ _children.remove( delChild );};
      
      
      std::list <Segment*> getChildren() { return _children;};
      std::list <Segment*> getParents()  { return _parents;};
      
      std::vector <IHit*> getHits()const {return _hits;};
      
      void addChild( Segment* child ){ _children.push_back(child); };
      void addParent( Segment* parent ){ _parents.push_back(parent); };
      
      unsigned getLayer()const { return _layer; };
      void setLayer( unsigned layer ) { _layer = layer; }; 
      
      std::vector<int>& getState() { return _state; }; 
      
      void raiseState() { if (_state.size() > 0) _state[0]++; };
      int getInnerState()const { return _state[0];}; 
      int getOuterState()const { return _state.back();}; 
      void resetState();
      
      void setSkippedLayers( unsigned skippedLayers ){ _state.resize( skippedLayers + 1 );}
      
      bool isActive() const { return _active;}
      void setActive( bool active ){ _active = active; }
      
      /** @return infos about the segment */
      std::string getInfo();
     
   private:
      
      std::list <Segment*> _children{}; 
      std::list <Segment*> _parents{};
      
      
      
      std::vector <IHit*> _hits{};
      
      std::vector<int> _state{};
      
      unsigned _layer{};
      bool _active{};
      
   };


} //end of KiTrack Namespace

#endif

