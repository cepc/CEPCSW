/* -*- Mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */
#ifndef NNClusterer_h
#define NNClusterer_h 1

#include <list>
#include <vector>

// #include "LCRTRelations.h"

/** Nearest neighbour type clusering for arbitrary types.
 *
 *  @author F.Gaede (DESY)
 *  @version $Id: NNClusterer.h 4052 2012-09-12 09:56:04Z gaede $
 */
namespace nnclu {

   /** fast check if integer is in a given range [Min,Max] */
  template <int Min, int Max >
  inline bool inRange( int i){   return ( (unsigned int) ( i - Min )  <= (unsigned int) ( Max - Min ) ); }

  /** fast check if integer is not in a given range [Min,Max] */
  template <int Min, int Max >
  inline bool notInRange( int i){   return ( (unsigned int) ( i - Min )  > (unsigned int) ( Max - Min ) ); }


  // forward declaration:
  template <class U>
  class Cluster ;


  /** Wrapper class for elements that are clustered, holding a pointer to the actual
   *  object "->first"  and a pointer to the cluster this obejct belongs to "->second".
   *
   *  @see Cluster
   *  @author F.Gaede (DESY)
   *  @version $Id: NNClusterer.h 4052 2012-09-12 09:56:04Z gaede $
   */
  template <class T>
  class Element : public  std::pair< T*, Cluster<T>* >{

    typedef T value_type ;
    typedef Cluster<T> cluster_type ;

  public:

    /** Default c'tor takes a pointer to the original element type object. The optional index can be used to
     *  code nearest neighbour bins, e.g. in z-coordinate to speed up the clustering process.
     */
    Element(T* element, int index0 = 0 ) : Index0( index0 )   {
      Pair::first = element ;
      Pair::second = 0 ;
    }

    /** C'tor that also takes a pointer to the cluster this element belongs to - in case seed elements/clusters are used.
     */
    Element(T* element ,  Cluster<T>* cl , int index0 = 0) : Index0( index0 ) {
      Pair::first =  element ;
      Pair::second = cl ;
    }

    /** Index that can be used to code nearest neighbour bins, e.g. in z-coordinate
     *  to speed up the clustering process.
     */
    int Index0 ;

  protected:
    typedef std::pair< T*, Cluster<T>* > Pair ;

    /** Don't allow default c'tor w/o element */
    Element() ;
  } ;


  /** Helper class that creates an Elements for an objects of type T.
   */
  template <class T>
  struct MakeElement{
    Element<T>  operator()( T* t) { return new Element<T>( t ) ;  }
  } ;


  /** Extension of std::vector that allows to take ownership of objects pointed to and
   *  delete these when going out of scope.
   */
  template <class T>
  class PtrVector : public std::vector<T*> {
    typedef std::vector<T*> vec ;
    bool _isOwner ;
  public:
    PtrVector() : _isOwner( false ) {}
    ~PtrVector() {
      if( _isOwner )
        for( typename vec::iterator i = vec::begin(),end = vec::end(); i != end ; delete *i++ ) ; //++i ) delete *i ;
    }
    void setOwner( bool val=true ) { _isOwner = val ; }
  };


  /** Extension of std::list that allows to take ownership of objects pointed to and
   *  delete these when going out of scope.
   */
  template <class T>
  class PtrList : public std::list<T*> {
    typedef std::list<T*> vec ;
    bool _isOwner ;
  public:
    PtrList() : _isOwner( false ) {}
    ~PtrList() {
      if( _isOwner )
        for( typename vec::iterator i = vec::begin(),end = vec::end(); i != end ; delete *i++ ) ; //++i ) delete *i ;
    }
    void setOwner( bool val=true ) { _isOwner = val ; }
  };



  /** Templated class for generic clusters  of Elements that are clustered with
   *  an NN-like clustering algorithm. Effectively this is just a list of elements.
   *
   *  @see Element
   *  @author F.Gaede (DESY)
   *  @version $Id: NNClusterer.h 4052 2012-09-12 09:56:04Z gaede $
   */
  template <class T >
  class Cluster : public std::list< Element<T> * >{ //, public lcrtrel::LCRTRelations {

  public :
    typedef Element<T> element_type ;
    typedef std::list< Element<T> * > base ;

    int ID ; //DEBUG

    Cluster() : ID(0) {}

    /** C'tor that takes the first element */
    Cluster( Element<T>* element)  {
      static int SID=0 ;  //DEBUG
      ID = SID++ ;      //DEBUG
      addElement( element ) ;
    }

    /** Add a element to this cluster - updates the element's pointer to cluster */
    void addElement( Element<T>* element ) {

      element->second = this ;
      base::push_back( element ) ;
    }

    // /** Remove all elements from the cluster and reset the cluster association, i.e. elements can be
    //  *  used for another clustering procedure.
    //  */
    // template <class Out>
    // void takeElements(Out result){
    //   typename Cluster<T>::iterator it = this->begin() ;
    //   while( it !=  this->end() ){
    //     (*it)->second = 0 ;
    //     result++ = *it ;
    //     it = this->erase(it) ;
    //   }
    // }


    /** Free all elements, ie. reset their cluster association - the elements are still in the list !.
     */
    void freeElements(){

      for( typename Cluster<T>::iterator it = this->begin(), end = this->end() ; it != end ; it++ ){
        (*it)->second = 0 ;
      }

      // typename Cluster<T>::iterator it = this->begin() ;
      // while( it !=  this->end() ){
      //   (*it)->second = 0 ;
      //   it = this->erase(it) ;
      // }
    }

    /** Merges all elements from the other cluster cl into this cluster */
    void mergeClusters( Cluster<T>* cl ) {

      for( typename Cluster<T>::iterator it = cl->begin(), end = cl->end() ; it != end ; it++ ){
        (*it)->second = this  ;
      }
      this->merge( *cl ) ;
    }

    /** D'tor frees all remaining elements that still belong to this cluster */
    ~Cluster()  {

      //  typename Cluster<T>::iterator it = this->begin() ;
      //  while( it !=  this->end()  )

      for( typename Cluster<T>::iterator it = this->begin() , end =  this->end()  ;  it != end ; ++it ){

        typename Cluster<T>::value_type h = *it ;

        if( h != 0 && h->second == this )
          h->second = 0 ;

      }
    }
  } ;


  template <class T>
  /** Main class for a nearest neighbour type clustering.
   *
   *  @author F.Gaede (DESY)
   *  @version $Id: NNClusterer.h 4052 2012-09-12 09:56:04Z gaede $
   */
  class NNClusterer{

  public:
    typedef T value_type ;
    typedef Cluster<T> cluster_type ;
    typedef Element<T> element_type ;
    typedef PtrVector< element_type >    element_vector ;
    typedef PtrVector< cluster_type >    cluster_vector ;
    typedef PtrList< element_type >      element_list ;
    typedef PtrList< cluster_type >      cluster_list ;

    /** Simple nearest neighbour (NN) clustering algorithm. Users have to provide an input iterator of
     *  Element objects and an output iterator for the clusters found. The predicate has to have
     *  a method with the following signature: bool operator()( const Element<T>*, const Element<T>*).
     *  All pairs of elements for which this method returns 'true' will be merged into one output cluster
     *  - all other pairs of elements will be in different clusters.
     */

    template <class In, class Out, class Pred >
    void cluster( In first, In last, Out result, Pred& pred , const unsigned minSize=1) {

      cluster_vector tmp ;
      tmp.reserve( 1024 ) ;

      while( first != last ) {

        for( In other = first+1 ;   other != last ; other ++ ) {

          if( pred( (*first) , (*other) ) ) {

            if( (*first)->second == 0 && (*other)->second == 0 ) {  // no cluster exists

              cluster_type* cl = new cluster_type( (*first) ) ;

              cl->addElement( (*other) ) ;

              tmp.push_back( cl ) ;

            }
            else if( (*first)->second != 0 && (*other)->second != 0 ) { // two clusters

              if(  (*first)->second != (*other)->second )  // don't call merge on identical clusters
                (*first)->second->mergeClusters( (*other)->second ) ;

            } else {  // one cluster exists

              if( (*first)->second != 0 ) {

                (*first)->second->addElement( (*other)  ) ;

              } else {

                (*other)->second->addElement( (*first)  ) ;
              }
            }

          }
        }
        ++first ;
      }

      // remove empty clusters
      for( typename cluster_vector::iterator i = tmp.begin(); i !=  tmp.end() ; i++ ){

        if( (*i)->size() > minSize-1 ) {

          result++ = *i ;
        }
        else {

          delete *i ;
        }
      }
    }


    /** Same as above - but requires the elements to be sorted in index0 (only compare neighbouring bins in index0). */

    template <class In, class Out, class Pred >
    void cluster_sorted( In first, In last, Out result, Pred& pred , const unsigned minSize=1) {


      cluster_vector tmp ;
      tmp.reserve( 1024 ) ;

      while( first != last ) {

        for( In other = first+1 ;   other != last ; other ++ ) {

          // if the elements are sorted we can skip the rest of the inner loop
          if( notInRange<-1,1>(   (*first)->Index0 - (*other)->Index0  )   )
            break ;

          if( pred( (*first) , (*other) ) ) {

            if( (*first)->second == 0 && (*other)->second == 0 ) {  // no cluster exists

              cluster_type* cl = new cluster_type( (*first) ) ;

              cl->addElement( (*other) ) ;

              tmp.push_back( cl ) ;

            }
            else if( (*first)->second != 0 && (*other)->second != 0 ) { // two clusters

              if(  (*first)->second != (*other)->second )  // don't call merge on identical clusters
                (*first)->second->mergeClusters( (*other)->second ) ;

            } else {  // one cluster exists

              if( (*first)->second != 0 ) {

                (*first)->second->addElement( (*other)  ) ;

              } else {

                (*other)->second->addElement( (*first)  ) ;
              }
            }

          }
        }
        ++first ;
      }

      // remove empty clusters
      for( typename cluster_vector::iterator i = tmp.begin(); i !=  tmp.end() ; i++ ){

        if( (*i)->size() > minSize-1 ) {

          result++ = *i ;
        }
        else {

          delete *i ;
        }
      }
    }

  };
  //-----------------------------------------------------------------------------------------------------------------------



  /**Splits a list into two based on a predicate. The new list will
   * hold all elements for which Pred is true which are in turn removed
   * from the original list.
   */
  template <class List, class Out, class Pred >
  void split_list( List& list, Out result, Pred pred) {

    typename List::iterator it = list.begin() ;

    while(  it != list.end() ){

      if( pred( *it ) ){

        result++ = *it ;
        it = list.erase( it ) ;

      } else
        ++it ;
    }
  }




  //-----------------------------------------------------------------------------------------------------------------------

} // namespace nnclu

#endif


