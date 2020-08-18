#include "KiTrack/Automaton.h"

#include <iostream>
//#include "marlin/VerbosityLevels.h"

using namespace KiTrack;

void Automaton::addSegment ( Segment* segment ){
  if ( segment->getLayer() >= _segments.size() ) { //in case this layer is not included so far
    _segments.resize( segment->getLayer() + 1 ); //resize the vector, so that the layer of the segment is now included
  }

  _segments[ segment->getLayer() ].push_back ( segment );

  _nConnections += segment->getChildren().size();
}

void Automaton::lengthenSegments(){

  // Info A: On skipped layers
  // ^^^^^^^^^^^^^^^^^^^^^^^^^^
  //
  // (read this only if you are interested on how the number of skipped layers is determined)
  //
  // The skipped layers are always between the innermost two hits of a segment.
  // Why? Because the connection between those two hits is the thing that differs
  // from its parent. Let's assume two 5-segments like those:
  //           /                                                               //
  //           \\                                                              //
  //           //                                                              //
  //           \\      _layer 2                                                //
  //            /      _layer 1                                                //
  // You see how they overlap with all of their hits except the inner one of the child
  // and the outer one of the parent. As the layer they are on equals the layer of the
  // innermost point, the child will have layer 1 and the parent layer 2.
  // (The choice of assigning the layer number of the innermost point and not the outermost
  // is arbitrary, it could be the other way as well.)
  //
  // Now suppose the child skips layer 1: ( in this example a kink in the segment means the layer is hit,
  // no kink means it is left out.)
  //
  //           /                                                               //
  //           \\                                                              //
  //           //                                                              //
  //           \\      _layer 2                                                //
  //            /      _layer 1                                                //
  //           /       _layer 0                                                //
  //
  // That means the child has now layer 0 and the parent layer 2.
  // This is no problem, the segment class has an outer and an inner state (simulating skipped layers)
  // instead of just an int. (to be more precise it has a vector containing the inner state and every layer left out)
  // Now when we want to make 6-segments (I know the numbers are high, but they help visualising), we would connect
  // parent and child to a new segment.
  //
  //
  //           /                                                               //
  //           \                                                               //
  //           /                                                               //
  //           \      _layer 2                                                 //
  //           /      _layer 1                                                 //
  //          /       _layer 0                                                 //
  //
  // So how many layers does this track skip: again 1 layer, so we need a state vector with 2 elements.
  // (one for layer 0, one for layer 1)
  //
  // So we don't care if there are any other skipped layers in the outer part of the segment, we only care about that
  // bit that won't overlap with parents.
  // So the easy recipe for the number of skipped layers after making a segment longer is:
  //   Compare the layer before ( 2 ) to the layer after ( 0 ). The skipped layers are the difference -1
  //   ( 2 - 0 - 1 = 1 --> segment->setSkippedLayers( 1 );
  
  //std::cout << "Combining the shorter segments to longer ones\n";
  
  //----------------------------------------------------------------------------------------------//
  //                                                                                              //
  // first: we create a new vector of a list of segments                                          //
  //   to have somewhere we can put the longer segments                                           //
  //                                                                                              //
  //----------------------------------------------------------------------------------------------//

  std::vector < std::list < Segment* > > longerSegments;
   
  if( _segments.size() > 0 ) longerSegments.resize ( _segments.size() -1 ); //This will have one layer less  

  //----------------------------------------------------------------------------------------------//
  //                                                                                              //
  // next: find all the longer segments and store them in the new vector[][]                      //
  //                                                                                              //
  //----------------------------------------------------------------------------------------------//
  
  unsigned nLongerSegments=0;
  unsigned nShorterSegments= _segments[0].size();
  
  for (unsigned layer = 1; layer < _segments.size(); layer++){ //over all layers where there still can be something below
    std::list<Segment*> segments = _segments[layer];
    for ( std::list<Segment*>::iterator iSeg=segments.begin(); iSeg != segments.end(); iSeg++){ //over all segments in this layer
      nShorterSegments++;
         
      Segment* parent = *iSeg;
      
      std::list <Segment*> children = parent->getChildren();
         
      for ( std::list<Segment*>::iterator iChild=children.begin(); iChild !=children.end(); iChild++){ //over all children of this parent
	Segment* child = *iChild;
            
	//Combine the parent and the child to form a new longer segment
        
	//take all the hits from the parent
	std::vector < IHit* > hits = parent->getHits();
        
	//and also add the inner hit from the child
	hits.insert( hits.begin(), child->getHits().at(0) );
        
	//make the new (longer) segment
	Segment* newSegment = new Segment ( hits );
	nLongerSegments++;
        
	//set the layer to the layer of the childsegment
	unsigned newLayer = child->getLayer();
	newSegment->setLayer ( newLayer );
        
	// Set the skipped layers.                  For an explanation see Info A above
	int skippedLayers = parent->getLayer() - child->getLayer() - 1;
	if( skippedLayers < 0 ) throw InvalidParameter( "skippedLayers can't be < 0!" );
	newSegment->setSkippedLayers( unsigned(skippedLayers) );      //
        /*
	std::cout << "Created longer segment: " << parent->getHits().size()
		  << "hits -->" << newSegment->getHits().size()
		  << " hits, layer = " << newLayer
		  << ", skipped layers = " << skippedLayers <<"\n";
	std::cout << "Combined: " << child->getInfo() << "<--with-->" << parent->getInfo() << "\n";
	*/
	// In a next step we want to again establish the conenctions between the longer segments (so we can do the 
	// Automaton and later combine them and then do it all again... ).
	// If we just created the Segments and dumped the old ones, we would have no idea what of the new, longer
	// Segments we can connect.
	// We could add some other container to store the possible connections of the longer Segments, but maybe
	// it's the easiest approach to use, what is already there: the shorter Segments.
	//
	// So when we combine two shorter segments, we store the new longer Segment as a parent or child.
	// Child, when the longer Segment goes on towards the inside, Parent if it continues on to the outside.
	// So the shorter Segments kind of act as joints, that hold the longer Segments together.
	//
	// Let's visulaize that, so that it makes more sense:
	// Let's have a look at 3 2-hit segments:
	//
	//          /       2-hit-Segment A
	//          \       2-hit-Segment B
	//          /       2-hit-Segment C
	//
	// Obviously we can make 2 3-hit segments out of this:
	//
	//          / -->   /       3-hit-Segment D
	//          \       \  \    .
	//          / -->      /    3-hit-Segment E
	// 
	// In the 2-hit-Segment B we store the 3-hit-Segments D and E as parent and child (while deleting A and C
	// as parent and child, because that is now not needed anymore )
	//
	// So when we want to connect the 3-hit-Segments, all we have to do is iterate over all 2-hit-Segments which
	// then only have 3-hit-Segments as parents and children.
	// When we come to Segment B, we see that D is a parent and E is a child, thus we connect them. Or to be more
	// precise, we connect them, if the criteria do say so.
	//
	// So, yes B acts like a joint connecting D and E
	// 
	// Erase the connection from the child to the parent segment and replace it with a link to the new
	// (longer) segment. ( and vice versa ) This way we can connect the longer segments easily after.
	child->deleteParent( parent );
	child->addParent ( newSegment );
	parent->deleteChild ( child );
	parent->addChild ( newSegment );
	// So now the new longer segment is a child of the old parent and a parent of the childsegment.
        
        
	// Save the new segment in the new vector[][]
	longerSegments[newLayer].push_back( newSegment );
      }
    }
  }

  //std::cout << " Made " << nLongerSegments << " longer segments from " << nShorterSegments << " shorter segments.\n";
  
  //----------------------------------------------------------------------------------------------//
  //                                                                                              //
  // Connect the new (longer) segments                                                            //
  //                                                                                              //
  //----------------------------------------------------------------------------------------------//
  
  //std::cout << "Next connecting the new longer segments\n";
     
  unsigned nConnections=0;
  unsigned nPossibleConnections=0;
  
  for ( unsigned layer = 1; layer + 1 < _segments.size(); layer++ ){ // over all layers (of course the first and the last ones are spared out because there is nothing more above or below
    std::list<Segment*> segments = _segments[layer];
    for ( std::list<Segment*>::iterator iSeg=segments.begin(); iSeg != segments.end(); iSeg++ ){ //over all (short) segments in this layer
      Segment* segment = *iSeg;
         
      std::list<Segment*> parents = segment->getParents();
      std::list<Segment*> children = segment->getChildren();
      
      for ( std::list<Segment*>::iterator iParent = parents.begin(); iParent != parents.end(); iParent++ ){ // over all parents of the segment
	Segment* parent = *iParent;
	for ( std::list<Segment*>::iterator iChild = children.begin(); iChild != children.end(); iChild++ ){ // over all children of the segment
	  Segment* child = *iChild;
	  
	  // Check if they are compatible
	  bool areCompatible = true;
	  ICriterion* theFailedCrit = NULL; 
          
	  //check all criteria (or at least until one returns false)
	  for ( unsigned iCrit = 0; iCrit < _criteria.size(); iCrit++ ){
	    if ( _criteria[iCrit]->areCompatible ( parent , child ) == false ){
	      areCompatible = false;
	      theFailedCrit = _criteria[iCrit];
	      break;
	    }
	  }
               
	  if ( areCompatible ){
	    //connect parent and child (i.e. connect the longer segments we previously created)
	    child->addParent( parent );
	    parent->addChild( child );
            
	    nConnections++;
            
	    //std::cout << "Connected: " << child->getInfo() << "<--with-->" << parent->getInfo() << "\n";
	  }
	  else{
	    //std::cout << "NOT Connected: " << child->getInfo() << "<--XXXX-->" << parent->getInfo() << "\n";
	    //if( theFailedCrit != NULL ) std::cout << "Failed first at criterion: " << theFailedCrit->getName() << "\n";
	  }
                         
	  nPossibleConnections++;
	}
      }
    }
  }
  _nConnections = nConnections;

  //std::cout << "Made " << nConnections << " of " << nPossibleConnections << " possible connections \n";
  
  //----------------------------------------------------------------------------------------------//
  //                                                                                              //
  //   Finally: replace the vector<list<segment*>> of the old segments with the new one           //
  //                                                                                              //
  //----------------------------------------------------------------------------------------------//
  
  //delete all old Segments:
  for( unsigned i=0; i<_segments.size(); i++){
    std::list<Segment*>& segments = _segments[i];
    for ( std::list<Segment*>::iterator iSeg=segments.begin(); iSeg != segments.end(); iSeg++ ){
      Segment* segment = *iSeg;
      delete segment;
    }
    segments.clear();
  }
   
  // And replace with the newer ones
  _segments = longerSegments;
}

void Automaton::doAutomaton(){

  bool hasChanged = true;
  int nIterations = -1;
  
  while ( hasChanged == true ){ //repeat this until no more changes happen (this should always be equal or smaller to the number of layers - 1
    hasChanged = false;
    nIterations++;
    
    for ( int layer = _segments.size()-1; layer >= 0; layer--){ //for all layers from outside in
      std::list <Segment*> segments = _segments[layer];
      for ( std::list<Segment*>::iterator iSeg=segments.begin(); iSeg!= segments.end(); iSeg++ ){ //for all segments in the layer
	Segment* parent= *iSeg;
            
	//Simulate skipped layers
	std::vector < int >& state = parent->getState();
	
	for ( int j= state.size()-1; j>=1; j--){
	  if ( state[j] == state[j-1] ){
	    state[j]++;
	    hasChanged = true; //something changed
	  }
	}
	
	if ( parent->isActive() ){
	  bool isActive = false; //whether the segment is active (i.e. still changing). This will be changed in the for loop, if it is active
	  
	  //Check if there is a neighbor
	  std::list <Segment*> children = parent->getChildren();
                         
	  for ( std::list<Segment*>::iterator iChild=children.begin(); iChild != children.end(); iChild++ ){// for all children
	    Segment* child = *iChild;

	    if ( child->getOuterState() == parent->getInnerState() ){  //Only if they have the same state
	      parent->raiseState(); //So it has a neighbor --> raise the state

	      hasChanged = true; //something changed
	      isActive = true;
	      
	      break; //It has a neighbor, we raised the state, so we need not check again in this iteration
	    }
	  }
               
	  parent->setActive( isActive );
	}
      }
    }
  }

  //std::cout << "Automaton performed using " << nIterations << " iterations.\n";
}

void Automaton::cleanBadStates(){

  unsigned nErasedSegments = 0;
  unsigned nKeptSegments = 0;

  for( unsigned layer=0; layer < _segments.size(); layer++ ){//for every layer
    std::list <Segment*> & segments = _segments[layer]; // We want to change things in the original list! Therefore the reference operator
    for( std::list<Segment*>::iterator iSeg= segments.begin(); iSeg != segments.end(); iSeg++ ){//over every segment
      Segment* segment = *iSeg;

      if( segment->getInnerState() == (int) layer ){ //the state is alright (equals the layer), this segment is good
	nKeptSegments++;
      }
      else { //state is wrong, delete the segment
	nErasedSegments++;

	//erase it from all its children
	std::list <Segment*> children = segment->getChildren();
	
	for (std::list<Segment*>::iterator iChild = children.begin(); iChild != children.end(); iChild++ ){
	  (*iChild)->deleteParent ( segment );
	  _nConnections--;
	}

	//erase it from all its parents
	std::list <Segment*> parents = segment->getParents();

	for (std::list<Segment*>::iterator iParent = parents.begin(); iParent!= parents.end(); iParent++){
	  (*iParent)->deleteChild ( segment );
	  _nConnections--;
	}

	//erase from the automaton
	delete *iSeg;
	iSeg = segments.erase( iSeg ); // erase the segment and update the iterator (updating is important!!!)
      }
    }
  }

  //std::cout << "Erased segments because of bad states= " << nErasedSegments << "\n";
  //std::cout << "Kept segments because of good states= " << nKeptSegments << "\n";
}

void Automaton::resetStates(){
  for ( unsigned layer = 0; layer < _segments.size(); layer++ ){ //over all layers
    std::list<Segment*> segments = _segments[layer];
    
    for ( std::list<Segment*>::iterator iSeg = segments.begin(); iSeg != segments.end(); iSeg++ ){ //over all segments in the layer
      (*iSeg)->resetState();
      (*iSeg)->setActive( true );
    }
  }
}

void Automaton::cleanBadConnections(){

  unsigned nConnectionsKept = 0;
  unsigned nConnectionsErased = 0;
  
  for ( int layer = _segments.size()-1 ; layer >= 1 ; layer-- ){ //over all layers from outside in. And there's no need to check layer 0, as it has no children.
    std::list<Segment*> segments = _segments[layer];
    for ( std::list<Segment*>::iterator iSeg = segments.begin(); iSeg!=segments.end(); iSeg++ ){ // over all segments in the layer
      Segment* parent = *iSeg;
      std::list < Segment* > children = parent->getChildren();
      
      for ( std::list<Segment*>::iterator iChild = children.begin(); iChild != children.end(); iChild++ ){ //over all children the segment has got
	Segment* child = *iChild;
	
	bool areCompatible = true; //whether segment and child are compatible
	
	//check all criteria (or at least until the first false pops up)
	for ( unsigned iCrit=0; iCrit < _criteria.size() ; iCrit++ ){
	  if ( _criteria[iCrit]->areCompatible( parent , child ) == false ){
	    areCompatible = false;
	    break; //no need to continue, now that we know, they're not compatible
	  }
	}
	
	if ( areCompatible == false ){ // they are not compatible --> erase the connection
	  nConnectionsErased++;
	  _nConnections--;
          
	  //erase the connection:
	  parent->deleteChild ( child );
	  child->deleteParent ( parent );
          
	  //A small note here: although we deleted a child from the vector, this doesn't mean we have to do iSeg--!
	  //Because we copied the value of segment->getChildren to the vector children. And this one doesn't change!
	}
	else{
	  nConnectionsKept++;
	}
      }
    }
  }
  
  //std::cout << "Erased bad connections= " << nConnectionsErased << "\n";
  //std::cout << "Kept good connections= " << nConnectionsKept << "\n";
}

std::vector < std::vector< IHit* > > Automaton::getTracksOfSegment ( Segment* segment, std::vector< IHit*> hits , unsigned minHits ){

  std::vector < std::vector< IHit* > > tracks; //the vector of the tracks to be returned

  std::vector <IHit*> segHits = segment->getHits(); // the hits of the segment

  //add the outer hit
  if ( segHits.back()->isVirtual() == false ) hits.push_back ( segHits.back() );  //Of course add only real hits to the track
  
  std::list <Segment*> children = segment->getChildren();

  if ( children.empty() ){ //No more children --> we are at the bottom --> start a new Track here
    //add the rest of the hits to the vector
    for ( int i = segHits.size()-2 ; i >= 0; i--){
      if ( segHits[i]->isVirtual() == false ) hits.push_back ( segHits[i] );
    }
          
    if ( hits.size() >= minHits ){
      //add this to the tracks
      tracks.push_back ( hits );
    }
  }
  else{// there are still children below --> so just take all their tracks and do it again
    for ( std::list<Segment*>::iterator iChild=children.begin(); iChild!= children.end(); iChild++){ //for all children
      std::vector < std::vector< IHit* > > newTracks = getTracksOfSegment( *iChild , hits );
      for (unsigned int j=0; j < newTracks.size(); j++){//for all the tracks of the child
	tracks.push_back ( newTracks[j] );
      }
    }
  }

  return tracks;
}

std::vector < std::vector< IHit* > > Automaton::getTracks( unsigned minHits ){

  std::vector < std::vector< IHit* > > tracks;
  std::vector <IHit*> emptyHitVec;

  for ( unsigned layer = 0 ; layer < _segments.size() ; layer++ ){ //over all layers
    std::list<Segment*> segments = _segments[layer];
    for ( std::list<Segment*>::iterator iSeg = segments.begin(); iSeg != segments.end(); iSeg++ ){ //over all segments
      Segment* segment = *iSeg;
      // by fucd: comment "if" in new ILC version of Automaton, why?
      if ( segment->getParents().empty() ){ // if it has no parents it is the end of a possible track
	// get the tracks from the segment
	std::vector < std::vector< IHit* > > newTracks = getTracksOfSegment( segment , emptyHitVec , minHits );
	
	// and add them to the vector of all tracks
	tracks.insert( tracks.end() , newTracks.begin() , newTracks.end() );
      }
    }
  }
  return tracks;
}

std::vector <const Segment*> Automaton::getSegments() const{
   
  std::vector <const Segment*> segments;
  
  for( unsigned layer=0; layer < _segments.size(); layer++ ){
    segments.insert( segments.end() , _segments[layer].begin() , _segments[layer].end() );
  }
     
  return segments; 
}

Automaton::~Automaton(){
  //delete the segments
  for( unsigned layer=0; layer < _segments.size(); layer++){ //over all layers
    std::list<Segment*> segments = _segments[layer];
    for( std::list<Segment*>::iterator iSeg = segments.begin(); iSeg!=segments.end(); iSeg++ ){ //over all segments
      delete *iSeg;
    }
  }
}







