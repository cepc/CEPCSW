#ifndef HopfieldNeuralNet_h
#define HopfieldNeuralNet_h
 
 
#include <vector>

#include "KiTrackExceptions.h"

namespace KiTrack{

   /**
    * Represents a Hopfield Neural Network
    * 
    * See <a href="../SubsetHopfieldNN.pdf">this</a> for detailed info.
    * 
    * Author: Robin Glattauer, HEPHY
    */
   class HopfieldNeuralNet {


      public:
         
         /**
         * @param G A matrix of the correlations between the neurons. 
         * True means two neurons are incompatible. False means, they are
         * compatible. (the diagonal elements are 0 by definition and any entry there will be ignored and set to 0) 
         * 
         * @param QI A vector containing the qualtity indicators of the neurons (i.e. their power to amplify or
         * weaken other neurons). Quality should be indicated by a value between 0 and 1. 1 being the highest quality.
         * 
         * @param states The states of the neurons. Should be between 0 and 1.
         * 
         * @param omega Controls the influence of the quality indicator on the  activation of the neuron. Needs to be
         * between 0 and 1. 0 means, no influence from the quality of the neurons -> system tends to biggest compatible 
         * set. 1 means highest influence from the quality of the neurons -> the highest quality neurons tend to win.
         */
         HopfieldNeuralNet( std::vector < std::vector <bool> > G , std::vector < double > QI , std::vector < double > states , double omega) ;
               
               
         /** Does one iteration of the neuronal network.
         * 
         * \f$ \vec{y} = W \times \vec{state} + \vec{w_0} \f$
         * 
         * \f$ \vec{state}_{new} = activationFunction(\vec{y}) \f$
         * 
         * @return Whether the Neural Network is considered as stable
         */
         bool doIteration();      
            
         /**
         * Sets the temperature of the Neural Network (The HNN is cooled down in every iteration)
         */
         void setT    (double T)    { _T = T;};
         
         /**
         * Sets the temperature at infinity. The temperature will converge to this
         * value after infinite iterations.
         */
         void setTInf (double TInf) {_TInf = TInf;};
         
         /**
         * Set the threshhold value below which the HNN is seen as "stable". As long as any Neuron changes its state 
         * by a value bigger than this, the HNN is not considered stable. When all Neurons change their states so little,
         * that none of the changes exceeds this threshold, then the HNN is stable.
         */
         void setLimitForStable (double limit) { _limitForStable = limit; };
         
         
         /** @return the vector of the states
         */
         std::vector <double> getStates(){ return _States; };
         

         
      protected:
         
         
            
         /** the matrix of the weights*/
         std::vector < std::vector <double> > _W{};
         
         /** states describing how active a neuron is*/
         std::vector < double > _States{};
         
         
         std::vector < double > _w0{};
         
         /** temperature */
         double _T{};
         
         /** temperature after infinite iterations */
         double _TInf{};

         /** indicates if the neuronal network is stable.
         * this is true when the change after one iteration 
         * of any neuron is not bigger than the value _limitForStable.
         */
         bool _isStable{};   

         /** The upper limit for change of a neuron, if it should be considered stabel.*/
         double _limitForStable{};
         
         /** Omega controls the influence of the quality indicator on the activation of the neuron.
         */
         double _omega{};

         /** the order of the neurons to be updated. So it should of course reach from 0 to the number of neurons -1.
         * (4 , 2, 0  1, 3) will for example mean: update first the neuron 4, then the neuron 2, then 0 and so on
         */
         std::vector <unsigned> _order{};
         
         
         /** Calculates the activation function
         * 
         * @param state the state
         * @param T the temperature
         * @return the activation function corresponding to the input values: g(x) = 1/2* (1 + tanh( x/T ))
         */
         double activationFunction ( double state , double T );


   };


}


#endif

