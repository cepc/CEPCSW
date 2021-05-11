#include "KiTrack/SegmentBuilder.h"

// ----- include for verbosity dependend logging ---------
//#include "marlin/VerbosityLevels.h"
#include <iostream>

using namespace KiTrack;

SegmentBuilder::SegmentBuilder(  std::map< int , std::vector< IHit* > > map_sector_hits ): 
  _map_sector_hits( map_sector_hits){}

Automaton SegmentBuilder::get1SegAutomaton(){
  /**********************************************************************************************/
  /*                Create and fill a map for the segments                                      */
  /**********************************************************************************************/
  std::map< int , std::vector< IHit* > >::iterator itSecHit; // Sec = sector , Hit = hits
  std::map< int , std::vector< Segment* > > map_sector_segments;
  std::map< int , std::vector< Segment* > > ::iterator itSecSeg; // Sec = sector , Seg = segments
        
  unsigned nCreatedSegments=0;
     
  for ( itSecHit = _map_sector_hits.begin(); itSecHit!=_map_sector_hits.end(); itSecHit++ ){ //over all sectors
    // All the hits in the sector
    int sector = itSecHit->first;
    std::vector <IHit*> hits = itSecHit->second;
    for ( unsigned int i=0; i < hits.size(); i++ ){ //over every hit in the sector
      // create a Segment
      Segment* segment = new Segment( hits[i] );
      segment->setLayer( hits[i]->getLayer() );
      
      // Store the segment in its map
      map_sector_segments[sector].push_back( segment );
      
      nCreatedSegments++;
    }
  }
     
  //std::cout << " Number of created 1-segments: " << nCreatedSegments <<"\n";
  
  /**********************************************************************************************/
  /*                Now check all 1-Segments and connect them to others                         */
  /*                Afterwards store them in an Automaton                                       */
  /**********************************************************************************************/
  
  unsigned nConnections=0;
  unsigned nStoredSegments = 0;
  
  Automaton automaton;
  
  for ( itSecSeg = map_sector_segments.begin(); itSecSeg != map_sector_segments.end(); itSecSeg++ ){ // over all sectors
    // All the segments with one certain code
    int sector = itSecSeg->first;
    std::vector <Segment*> segments = itSecSeg->second;
          
    // Now find out, what the allowed codes to connect to are:
    std::set <int> targetSectors;
    
    for ( unsigned i=0; i < _sectorConnectors.size(); i++ ){ // over all IHitConnectors
      // get the allowed targets
      std::set <int> newTargetSectors = _sectorConnectors[i]->getTargetSectors( sector );
      
      //insert them into our set
      targetSectors.insert( newTargetSectors.begin() , newTargetSectors.end() );
    }
          
    for ( unsigned int i=0; i< segments.size(); i++ ){ //over all segments within the sector
      Segment* parent = segments[i]; 
      
      for ( std::set<int>::iterator itTarg = targetSectors.begin(); itTarg!=targetSectors.end(); itTarg++ ){ // over all target codes
	int targetSector = *itTarg;
	std::vector <Segment*> targetSegments = map_sector_segments[ targetSector ];
                    
	for ( unsigned int j=0; j < targetSegments.size(); j++ ){ // over all segments in the target sector
	  Segment* child = targetSegments[j];
	  bool areCompatible = true;
	  ICriterion* theFailedCrit = NULL; 
                         
	  for (unsigned int iCrit = 0; iCrit < _criteria.size(); iCrit++){
	    if ( _criteria[iCrit]->areCompatible( parent , child ) == false ){
	      areCompatible = false;
	      theFailedCrit = _criteria[iCrit];
	      break;
	    }
	  }
                         
	  if ( areCompatible ){ //the connection was successful 
	    parent->addChild( child );
	    child->addParent( parent );
            
	    nConnections++;              
	    //std::cout << "Connected: " << child->getInfo() << "<--with-->" << parent->getInfo() << "\n"; 
	  } 
	  else{
	    //std::cout << "NOT Connected: " << child->getInfo() << "<--XXXX-->" << parent->getInfo() << "\n";
	    //if( theFailedCrit != NULL ) std::cout << "Failed first at criterion: " << theFailedCrit->getName() << "\n";
	  }
	}
      }
         
      // Store the segment in the automaton
      automaton.addSegment( parent );
      nStoredSegments++;
    }      
  }
      
  //std::cout << "Number of connections made " << nConnections <<"\n";
  //std::cout << "Number of 1-segments, that got stored in the automaton: " << nStoredSegments <<"\n";
    
  return automaton;
}







