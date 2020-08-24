#include "KiTrack/HopfieldNeuralNet.h"

#include <cmath>
#include <iostream>
#include <algorithm>
#include <sstream>
#include <random>


using namespace KiTrack;

HopfieldNeuralNet::HopfieldNeuralNet( std::vector < std::vector <bool> > G , std::vector < double > QI , std::vector < double > states , double omega) {

   unsigned int nNeurons = G.size();

   
   
   //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   // Check the validity of the input parameters
   
   
   std::stringstream s;
   s << "HopfieldNeuralNet: ";
   
   // Is G a square matrix?
   for( unsigned i=0; i< nNeurons; i++ ){
      
      if( G[i].size() != nNeurons ){
         
         s << "G must be a square matrix!!! G[" << i << "].size() == " << G[i].size() << " != G.size() (" << nNeurons << ")\n";
         throw InvalidParameter( s.str() );
         
      }
      
   }
   
   // Does the QI vector have the right size?
   if( QI.size() != nNeurons ){
      
      s << "The QI vector must have the same size as G! QI.size() == " << QI.size() << " != G.size() (" << nNeurons << ")\n";
      throw InvalidParameter( s.str() );
      
   }
   
   // Are all the Quality Indicators in the range 0 - 1
   for( unsigned i=0; i< nNeurons; i++ ){
      
      if( ( QI[i] < 0. ) || ( QI[i] > 1. ) ){
         
         s << "The QI must be between 0 and 1, QI[" << i << "] == " << QI[i] << " is not a valid value!\n";
         
      }
      
   }
   
   // Does the states vector have the right size?
   if( states.size() != nNeurons ){
      
      s << "The vector of the states must have the same size as G! states.size() == " << states.size() << " != G.size() (" << nNeurons << ")\n";
      throw InvalidParameter( s.str() );
      
   }
   
   // Are all the states in the range 0 - 1
   for( unsigned i=0; i< nNeurons; i++ ){
      
      if( ( states[i] < 0. ) || ( states[i] > 1. ) ){
         
         s << "The states must be between 0 and 1, states[" << i << "] == " << states[i] << " is not a valid value!\n";
         
      }
      
   }
   
   // Is omega in the range from  0 to 1
   if( ( omega < 0. ) || ( omega > 1. ) ){
      
      s << "Omega must be in the range from 0 to 1, omega == " << omega << " is not a valid value!\n";
      
   }
   // End of checking the parameters
   //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   
   _omega = omega;
   _States = states;
   
   // resize the vectors.
   _States.resize( nNeurons );
   _w0.resize( nNeurons );
   _W.resize( nNeurons );
   _order.resize( nNeurons);
   
   for ( unsigned int i =0; i < nNeurons; i++){
      
      // resize the vectors of the matrix W
      _W[i].resize( nNeurons );
      
      
      // initialise the order vector
      _order[i]=i;                      //the order now is 0,1,2,3... (will be changed to a random sequence in the iteration)
      
   }
   
   
   //calculate _w0
   for (unsigned int i=0; i < QI.size(); i++) _w0[i] = omega * QI[i];
   
   
   
   // Build the W matrix. (the matrix of the influences including their force)
   
   double comp = 1;
   if (nNeurons > 0 ) comp = (1. - omega) / double (nNeurons);
   

   for (unsigned int i=0; i< nNeurons ; i++){ 

      for (unsigned int j=0; j< nNeurons ; j++){
       
         if (i == j) _W[i][j] = 0.; //diagonal elements are 0 --> whatever the matrix G says here is ignored.
            
         else{  
            
            if ( G[i][j] == 1 ) _W[i][j] = -1;   //Neurons are incompatible
            
            else _W[i][j] =  comp;   //Neurons are compatible
            
         }
         
      }
      
   }
   
  
   


   _T = 0;
   _TInf = 0;
   
   _isStable = false;
   _limitForStable = 0.01;



}



double HopfieldNeuralNet::activationFunction ( double state , double T ){
   
   
   double y = 1;
   
   if (T > 0) y = 0.5 * ( 1 + tanh( state / T ) ); //if T==0  tanh( infinity ) gives 1.


   return y;
   
   
}




bool HopfieldNeuralNet::doIteration(){
   
   _isStable = true;
   
   // initialize the random generator
   std::random_device rng;
   std::mt19937 urng(rng());

   shuffle ( _order.begin() , _order.end() , urng ); //shuffle the order
   
   for (unsigned int i=0; i<_States.size() ; i++){ //for all entries of the vector
      
      unsigned iNeuron = _order[i];
      
      
      double y;
      
      y = _w0[iNeuron];
      
      //matrix vector multiplication (or one line of it to be precise)  
      for (unsigned int j=0; j< _W[iNeuron].size(); j++){ 
       
         y  += _W[iNeuron][j] * _States[j]; 
         
      }
      
      y = activationFunction ( y , _T );
      
      // check if the change was big enough that the Network is not stable
      if ( fabs( _States[iNeuron] - y ) > _limitForStable ) _isStable = false;
      
      // update the state
      _States[iNeuron] = y;
      
   }
   
   
   
   // after the iteration, we need to calculate a new tempereatur _T
   
   _T = 0.5* ( _T + _TInf);
   
   
   return _isStable;
   
}






