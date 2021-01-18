#ifndef ConfigFlags_h
#define ConfigFlags_h

#include <string>
#include <iostream>
#include <exception>
#include <map>

namespace MarlinTrk {
  
  class ConfigFlags ;
  
  inline std::ostream& operator<<(std::ostream& os, const ConfigFlags& cf) ;
  
  class ConfigFlags{
    
    friend  std::ostream& operator<<(std::ostream& os, const ConfigFlags& flags) ;
    
    typedef std::pair<std::string, bool> Flag ;
    typedef std::map< unsigned, Flag > Map ;
    
    
  public:
    
    /** Helper class that holds a  number of boolean properties for configuration.
     *  The property keys  (type: unsigned) have to be defined in the class using 
     *  this class.
     */
    ConfigFlags() {}
    
    ~ConfigFlags(){}
    
    
    void registerOption( unsigned key, const std::string& name, bool defaultValue=false ){
      
      _map[ key ] = std::make_pair( name , defaultValue ) ;
    }
    
    bool option(unsigned key) const { 
      
      Map::const_iterator it = _map.find( key ) ;
      
      if( it == _map.end() )
        return false ;
      
      return it->second.second  ;
    }
    
    bool operator[](unsigned key) const { 
      return option( key ) ;
    }
    
    void setOption(unsigned key , bool val) { 
      
      Map::iterator it = _map.find( key ) ;
      
      if( it !=_map.end() )
        it->second.second  = val ;
    }
    
    
    std::string& optionName(unsigned key) {
      
      static std::string empty("UNKNOWN") ;
      
      Map::iterator it = _map.find( key ) ;
      
      if( it == _map.end() ) 
        return empty ; 
      
      return it->second.first ;
    }
    
  protected:
    Map _map ;
    
  };
  
  
  inline std::ostream& operator<<(std::ostream& os, const MarlinTrk::ConfigFlags& cf)  {
    
    for( ConfigFlags::Map::const_iterator it = cf._map.begin(); it != cf._map.end()  ; ++it){
      
      os << "  option: " << it->second.first <<  "\t: " << it->second.second << std::endl ; 
    }
    
    return os ;
  } 
  
} // namespace


#endif

