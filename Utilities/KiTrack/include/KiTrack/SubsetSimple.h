#ifndef SubsetSimple_h
#define SubsetSimple_h

#include <algorithm>

//#include "marlin/VerbosityLevels.h"

#include "KiTrack/Subset.h"

namespace KiTrack {
   
   
   template <class T, class QI >
   struct SorterQI
   {
      SorterQI( QI getQI ): _getQI( getQI ){}
      
      bool operator()( T a, T b ){ return ( _getQI( a ) > _getQI( b ) ); }
      
      QI _getQI;
   };
   
 
   /** A class to get the best subset by taking the element with the highest quality
    * indicator and throwing away all incompatible ones. 
    * On the remaining elements repeat the procedure. 
    * 
    */
   template <class T >
   class SubsetSimple : public Subset<T> {
      
   public:
      
      /** Calculates the best set using a very simple method. After this the results can be accessed with 
       * the getAccepted() and getRejected() methods (provided by the baseclass Subset).
       * 
       * This is a templated class to allow the reuse with objects from different.
       * The goal is this: if you have a set of things that are somehow compatible or incompatible 
       * with each other and you want a subset that only contains compatible ones and has a high quality, to find
       * this subset.
       * 
       * For example: you plan a concert where a lot of different artists shall perform together. You have 20 that
       * you are interested in asking. BUT: not all of them get along together very well. So artist 1 might be incompatible
       * with artist 2 and 7 for some arbitrary reasons. And artist 2 is incompatible with artist 1,3 and 13 and so on.
       * In order not to completely mess up the concert you can only invite artists who are entirely compatible with each other.
       * AND: not all artists are equal: some are really famous and some are pretty mediocre.
       * So you want to find a set of artists with the highest possible quality which is completely compatible.
       * 
       * This is done by this class. 
       * 
       * The algorithm used for it is this: take the element with the highest quality and accept it. Reject all
       * other elements that are incompatible with it. From the remaining elements redo the same with the next best
       * element.
       * In our artist example: take the biggest star and throw out all he doesn't get allong with. From the remaining
       * artists (those who do get along with the biggest star) pick next famous one. Once again kick out all that
       * don't get along with him and so on.
       * 
       * In order to work it needs to know two things: how to calculate the quality of an element and how to determine
       * if two elements are compatible. These are of course things the user has to provide.
       * For both a functor object is needed.
       * 
       * Here is how it could look like in the artist example:
       * 
\verbatim

class AristQI{
   
public:
   
   // returns the quality of an artist (a value between 0 and 1, 0 being bad and 1 being fantastic)  
   inline double operator()( Artist artist ){
      
      return artist.numberOfFans()/nPeopleOnEarth;    // a number between 0 and 1
      
   }
   
};

class ArtistCompatibility{
   
public:
   
   // returns whether two artists are compatible
   inline bool operator()( Artist artistA, Artist artistB ){
      
      if( artistA.hates( artistB ) ) return false;
      else return true; 
      
   }
   
};


//Somewhere within the program:

ArtistCompatibility comp;
ArtistQI qi;

SubsetHopfieldNN< Artist > subset;
subset.add( vecOfArtists );                                 
subset.calculateBestSet( comp, qi );

std::vector< Artist > artistsToInvite = subset.getAccepted();
std::vector< Artist > artistsToStayAtHome = subset.getRejected();


\endverbatim
       * 
       * 
       * 
       * @param areCompatible a functor of type bool( T, T ) that should tell whether two elements are compatible
       * or not.
       * 
       * @param getQI a functor of type double( T ) that returns the quality of an element and should range between 0 and 1.
       */
      template< class GetQI, class AreCompatible >
      void calculateBestSet( AreCompatible areCompatible, GetQI getQI );
   };
   
   template<class T> template<class GetQI, class AreCompatible> 
   void SubsetSimple<T>::calculateBestSet( AreCompatible areCompatible, GetQI getQI ){
      
      
      unsigned nAccepted=0;
      unsigned nRejected=0;
      
      
      
      std::vector< T > elements = this->_elements;
      
      // sort the vector from big QI to small QI
      SorterQI< T, GetQI > sorterQI( getQI );
      sort( elements.begin(), elements.end() , sorterQI );
      
      //std::cout << "DEBUG: SubsetSimple::calculateBestSet The elements and their QIs sorted:" << std::endl;
      for( unsigned i=0; i < elements.size(); i++ ){
	double qi = getQI( elements[i] );
	//std::cout << "DEBUG: SubsetSimple::calculateBestSet " << elements[i] << "\t" << qi << std::endl;
      }
      
      /*The idea here is: the first track is (now that we sorted) the one with the highest
       * QI. Check all other tracks if they are complatible with it.
       * If one of them is incompatible throw it away. Once done with that, store the
       * first track as accepted (and delete it from the vector) and do it once more for
       * the new highest QI track.
       * Repeat until all tracks are stored.
       */   
      while( elements.size() > 0 ){
         
         
         // check all elements with smaller QI if they are compatible
         for( unsigned i=1; i < elements.size(); i++){
            
            
            if( areCompatible( elements[0] , elements[i] ) == false ){ // the track is incompatible
               
               // reject it
               nRejected++;
               this->_rejectedElements.push_back( elements[i] );
               
	       //std::cout << "DEBUG: SubsetSimple::calculateBestSet " << "reject " << elements[i] << std::endl;
               // and delete it from the elements we are looking at
               elements.erase( elements.begin() + i );
               i--;
               
            }
            
         }
         
         // store the first element
	 //std::cout << "DEBUG: SubsetSimple::calculateBestSet " << "accept " << elements[0] << std::endl;
         nAccepted++;
         this->_acceptedElements.push_back( elements[0] );
         
         // delete it from the current tracks
         elements.erase( elements.begin() );
         
      }
            
      //std::cout << "DEBUG: SubsetSimple::calculateBestSet " << "SubsetSimple accepted " << nAccepted 
      //		<< " elements and rejected " << nRejected << " elements of all in all " 
      //		<< nAccepted + nRejected << " elements." << std::endl;
      
   }
} // end namespace
#endif
