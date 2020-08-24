#ifndef ICriterion_h
#define ICriterion_h

#include <vector>
#include <map>
#include <string>

#include "KiTrack/Segment.h"
#include "KiTrack/KiTrackExceptions.h"



namespace KiTrack{


   /** An Interface for Criteria.
    * 
    * A Criterion is a class, that is able to take two Segments and check whether they are compatible or not.
    */   
   class ICriterion{


   public: 
      
      /** @return If the two Segments are compatible with each other, i.e. could be combined to a longer Segment or a 
       * track.
       */
      virtual bool areCompatible( Segment* parent , Segment* child ) = 0;
      
      
      /** @return A map, where the calculated values are stored. The keys are the names of the values.
       */
      std::map < std::string , float > getMapOfValues() {return _map_name_value; };
      
      
      virtual ~ICriterion(){};
      
      /** Sets, whether the calculated values shall be saved in a map
       */
      void setSaveValues( bool saveValues ){ _saveValues = saveValues;}
      
      
      /** @return the name of the criterion */
      std::string getName(){return _name;}
      
      /** @return the type of the criterion */
      std::string getType(){return _type;}
     
   protected:
      
      
      std::map < std::string , float > _map_name_value{};
      
      bool _saveValues{};
      
      std::string _name{};
      std::string _type{};
      
   };
   
}


#endif




