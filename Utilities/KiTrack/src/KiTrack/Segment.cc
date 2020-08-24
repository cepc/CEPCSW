#include "KiTrack/Segment.h"

#include <sstream>

using namespace KiTrack;

Segment::Segment( std::vector <IHit*> hits){ 

   _hits = hits; 
   
   _state.push_back(0); 
   
   _children.clear(); 
   _parents.clear();

   _active = true;
   
   _layer=0;
}



Segment::Segment( IHit* hit){ 
   
   _hits.push_back( hit) ;
   _state.push_back(0); 
   _children.clear(); 
   _parents.clear();
   
   _active = true;
   
   _layer=0;
}



void Segment::resetState(){
   
   
   for ( unsigned i = 0; i <_state.size(); i++){
      
      _state[i] = 0;
      
   }
   
}

std::string Segment::getInfo(){
   
 
   std::stringstream info;
   
   for( unsigned i=0; i<_hits.size(); i++ ) info << _hits[i]->getPositionInfo();
   
   info << "[";
   for( unsigned i=0; i+1<_state.size(); i++ ) info << _state[i] << ",";
   info << _state.back() << "]";
   
   return info.str();
   
   
}



   

