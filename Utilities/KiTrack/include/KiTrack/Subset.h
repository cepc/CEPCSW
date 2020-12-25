#ifndef Subset_h
#define Subset_h


#include "KiTrack/ITrack.h"


namespace KiTrack {
 
   /** A base class for subsets.
    * 
    * The idea behind the subset classes is, that one has a collection of elements, but not all elements are compatible
    * with each other. One wants to find a subset of these elements, where all are compatible. Also the subset should have
    * a high quality. 
    * Details about this can be found in the derived classes. For technical reasons (templates and virtual methods) 
    * this class only covers the basics.
    */
   template <class T>
   class Subset{
      
   public:
      
      /** Adds an element
       */
      void add( T newElement ){ _elements.push_back( newElement ); };
      
      /** Adds a vector of elements
       */
      void add( std::vector<T> newElements ){ _elements.insert( _elements.end() , newElements.begin() , newElements.end() ); };
      
      
      /** @return the subset of the best tracks, i.e. all the tracks that were accepted
       */
      std::vector< T > getAccepted(){ return _acceptedElements;} ;
      
      /** @return the tracks that got rejected (and are therefore not in the best subset)
       */
      std::vector< T > getRejected(){ return _rejectedElements;} ;
      
      
   
   protected:
      
      std::vector< T > _elements{};
      std::vector< T > _acceptedElements{};
      std::vector< T > _rejectedElements{};
      
      
   };
   
   
   
   
}


#endif

