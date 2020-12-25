#ifndef ITrack_h
#define ITrack_h

#include <vector>
#include "KiTrack/IHit.h"

namespace KiTrack{

   
   /** Abstract Base Class for tracks.
    */
   class ITrack{
      
      
   public:
      
      /** @return a vector of the hits of the track
       */
      virtual std::vector< IHit*> getHits() const = 0;
      
      /** @return an indicator for the quality of the track. Usually between 0 and 1 */
      virtual double getQI() const = 0; 
      
      /** Fits the track. */
      virtual void fit() = 0;
      
      /** @return the degrees of freedom of the fit. This must only be called after fit() has been used */
      virtual double getNdf() const = 0;
      
      /** @return the chi squared value of the fit. This must only be called after fit() has been used */
      virtual double getChi2() const = 0;
      
      /** @return the chi squared probability of the fit. This must only be called after fit() has been used */
      virtual double getChi2Prob() const = 0;
      
      virtual ~ITrack(){}
      
   };




}

#endif


