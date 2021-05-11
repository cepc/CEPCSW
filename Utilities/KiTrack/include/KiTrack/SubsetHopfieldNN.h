#ifndef SubsetHopfieldNN_h
#define SubsetHopfieldNN_h

#include <CLHEP/Random/RandFlat.h>

#include "KiTrack/Subset.h"
#include "KiTrack/HopfieldNeuralNet.h"



namespace KiTrack {
 
   /** A class to get the best subset with help of a Hopfield Neural Network
    * 
    */
   template< class T >
   class SubsetHopfieldNN : public Subset<T>{
      
   public:
      
      /** Calculates the best set using a Hopfield Neural Network 
       * (see <a href="../SubsetHopfieldNN.pdf">this</a> for more info) . 
       * 
       * After this the results can be accessed with 
       * the getAccepted() and getRejected() methods (provided by the baseclass Subset).
       *
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
       * This is done by this class. The algorithm used for it is a Hopfield Neural Network (HNN).
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
       * @param areCompatible a functor of type bool( T, T ) that should tell whether two elements are compatible
       * or not.
       * 
       * @param getQI a functor of type double( T ) that returns the quality of an element and should range between 0 and 1.
       */
      template< class GetQI, class AreCompatible >
      void calculateBestSet( AreCompatible areCompatible, GetQI getQI );
      
      
      SubsetHopfieldNN(){ 
       
         _TStart = 2.1;
         _TInf = 0.1;
         _omega = 0.75;
         _limitForStable = 0.01;
         _initStateMin = 0.;
         _initStateMax = 0.1;
         _activationThreshold = 0.5;
         
      }
      
      void setTStart( double tStart ){ _TStart = tStart; }
      void setTInf( double tInf ){ _TInf = tInf; }
      void setOmega( double omega ){ _omega = omega; }
      void setLimitForStable( double limitForStable ){ _limitForStable = limitForStable; }
      void setInitStateMin( double initStateMin ){ _initStateMin = initStateMin; }
      void setInitStateMax( double initStateMax ){ _initStateMax = initStateMax; }
      void setActivationThreshold( double activationThreshold ){ _activationThreshold = activationThreshold; }
      
      double getTStart(){ return _TStart; }
      double getTInf(){ return _TInf; }
      double getOmega(){ return _omega; }
      double getLimitForStable(){ return _limitForStable; }
      double getInitStateMin(){ return _initStateMin; }
      double getInitStateMax(){ return _initStateMax; }
      double getActivationThreshold(){ return _activationThreshold; }
      
   protected:
      
      
      double _TStart{};
      double _TInf{};
      double _omega{};
      double _limitForStable{};
      double _initStateMin{};
      double _initStateMax{};
      double _activationThreshold{};
      
   };
   
   
   template< class T > template< class GetQI, class AreCompatible >
   void SubsetHopfieldNN<T>::calculateBestSet( AreCompatible areCompatible, GetQI getQI ){
      
      
      unsigned nAccepted=0;
      unsigned nRejected=0;
      unsigned nCompWithAll=0;
      unsigned nIncompatible=0;
      
      
      std::vector< T > elements = this->_elements; //this pointer is needed here, because of the template!
      
      unsigned nElements = elements.size();
      
      // the information for the Hopfield Neural Network:
      
      std::vector < std::vector <bool> > G; // a matrix telling, if two neurons (elements) are compatible
      G.resize( nElements );
      for (unsigned i=0; i<nElements; i++) G[i].resize( elements.size() );
      
      std::vector < double > QI ; // the quality indicators of the neurons (elements)
      QI.resize( nElements );
      
      std::vector < double > states; // the initial state to start from.
      states.resize( nElements );
      
      
      
      
      /**********************************************************************************************/
      /*                1. Find out which elements are compatible and get the QIs                     */
      /**********************************************************************************************/
      
      
      for ( unsigned i=0; i < nElements ; i++){ //over all elements
         
         
         T elementA = elements[i]; //the track we want to look at.
         
         // Get the quality
         QI[i] = getQI( elementA );
         
         //streamlog_out(DEBUG3) << "QI of element " << i << " = " << QI[i] << "\n";
         
         
         // Set an initial state
         states[i] = CLHEP::RandFlat::shoot ( _initStateMin , _initStateMax ); //random ( uniformly ) values from initStateMin to initStateMax
         
         
         // Fill the states in the G matrix. (whether two elements are compatible or not
         for ( unsigned j=i+1; j < nElements ; j++ ){ // over all elements that come after the current one (the elements before get filled automatically because of symmetry)
            
            T elementB = elements[j]; // the track we check if it is in conflict with trackA
            
            if ( areCompatible( elementA , elementB ) ){ 
               
               G[i][j] = 0;
               G[j][i] = 0;
               
            }
            else{
               
               G[i][j] = 1;
               G[j][i] = 1;            
               
            }
            
         }
         
      }
      
      // output of the G matrix:
      if( !G.empty() ){
         
	//streamlog_out(DEBUG2) << "G:\n";
         
         
         for ( unsigned i=0; i < G.size(); i++ ){
            
            
            for ( unsigned j=0; j < G[i].size(); j++ ){
               
	      //streamlog_out(DEBUG2) << G[i][j] << "  ";
               
            }
            
            //streamlog_out(DEBUG2) << "\n";
            
         }
        
      }
      // output, where one sees, what elements are  incompatible with what others:
      if( !G.empty() ){
         
	//streamlog_out(DEBUG2) << "Incompatible ones:\n";
         
         
         for ( unsigned i=0; i < G.size(); i++ ){
            
	   //streamlog_out(DEBUG2) << "Element " << i << ":\t";
            
            for ( unsigned j=0; j < G[i].size(); j++ ){
               
	      //if( G[i][j] ) streamlog_out(DEBUG2) << j << ", ";
               
            }
            
            //streamlog_out(DEBUG2) << "\n";
            
         }
         
      }
      
      
      /**********************************************************************************************/
      /*                2. Save elements, that are compatible with all others                         */
      /**********************************************************************************************/
      
      for( unsigned i=0; i < elements.size(); i++ ){
         
         
         bool isCompatibleWithAll = true;
         
         //check if this track is compatible with all others
         for( unsigned j=0; j < elements.size(); j++){
            
            //G[i][j] == 0 i and j are compatible, if G[i][j] == 1 i and j are incompatible
            if( (i != j) && G[i][j] ){
               
               isCompatibleWithAll = false;
               break;
               
            }
            
         }
         
         
         if ( isCompatibleWithAll ){ //if it is compatible with all others, we don't need the Hopfield Neural Net, we can just save it
            
            
            //add the track to the good ones
            this->_acceptedElements.push_back( elements[i] );
            nCompWithAll++;
            
            
            //And now erase it from the ones we will still check:
            elements.erase( elements.begin() + i );
            states.erase( states.begin() + i );
            QI.erase( QI.begin() + i );
            
            for( unsigned j=0; j<G.size(); j++ ) G[j].erase( G[j].begin() + i );
                     G.erase( G.begin() + i );
            
            i--;
            
         }
         else{
            
            nIncompatible++;
            
         }
         
      }
      
      
      //streamlog_out( DEBUG3 ) << nCompWithAll << " elements are compatible with all others, " << nIncompatible
      //<< " elements are interfering and will be checked for the best subset\n";
      
      
      
      /**********************************************************************************************/
      /*                3. Let the Neural Network perform to find the best subset                   */
      /**********************************************************************************************/  
      
      if( !elements.empty() ){
         
         HopfieldNeuralNet net( G , QI , states , _omega);
         
         net.setT ( _TStart );
         net.setTInf( _TInf );
         net.setLimitForStable( _limitForStable );
         
         unsigned nIterations=1;
         
         //streamlog_out(DEBUG1) << "states: ( ";
         //for ( unsigned int i=0; i< states.size(); i++) streamlog_out(DEBUG1) << states[i] << " "; 
         //streamlog_out(DEBUG1) << ")\n";
         
         while ( !net.doIteration() ){ // while the Neural Net is not (yet) stable
            
            nIterations++;
            
            std::vector <double> newStates = net.getStates();
            
            //streamlog_out(DEBUG1) << "states: ( ";      
            
            //for ( unsigned int i=0; i< newStates.size(); i++) streamlog_out(DEBUG1) << newStates[i] << " "; 
            
            //streamlog_out(DEBUG1) << ")\n";
            
         }
         
         
         
         //streamlog_out( DEBUG3 ) << "Hopfield Neural Network is stable after " << nIterations << " iterations.\n";
         
         
         
         
         /**********************************************************************************************/
         /*                4. Now just sort the elements into accepted and rejected ones                 */
         /**********************************************************************************************/  
         
         
         states = net.getStates();
         
         
         
         for ( unsigned i=0; i < states.size(); i++ ){
            
            
            if ( states[i] >= _activationThreshold ){
               
               this->_acceptedElements.push_back( elements[i] );
               nAccepted++;
               
            }
            else{
               
               this->_rejectedElements.push_back( elements[i] );
               nRejected++;
               
            }
            
         }
         
      }
      
      
      //streamlog_out( DEBUG3 ) << "Hopfield Neural Network accepted " << nAccepted 
      //<< " elements and rejected " << nRejected << " elements of all in all " 
      //<< nAccepted + nRejected << "incomaptible elements.\n";
      
      //streamlog_out( DEBUG3 )   << "So in sum " << nAccepted + nCompWithAll
      //<< " elements survived and " << nRejected << " elements got rejected.\n";
      
      
   }
   
   
}


#endif

