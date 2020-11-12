#ifndef Other_Helpers_hh
#define Other_Helpers_hh 1

#include <iostream>
#include <map>
#include <stdexcept>

namespace CEPC {
  typedef enum {
    kCenter                     = 0,
    kCenterSide                 = 1,
    kWaist                      = 2,
    kFatWaist                   = 3,
    kCrotch                     = 4,
    kCrotchAsymUp               = 5,
    kCrotchAsymDn               = 6,
    kLegs                       = 7,
    kFlareLegUp                 = 8,
    kFlareLegDn                 = 9
  } ECrossType;
  
  inline ECrossType getCrossType( std::string const & type) {
    std::map< std::string, CEPC::ECrossType > CrossTypes;
    CrossTypes["Center"]                = CEPC::kCenter;
    CrossTypes["CenterSide"]            = CEPC::kCenterSide;
    CrossTypes["Waist"]                 = CEPC::kWaist;
    CrossTypes["FatWaist"]              = CEPC::kFatWaist;
    CrossTypes["Crotch"]                = CEPC::kCrotch;
    CrossTypes["CrotchAsymUp"]          = CEPC::kCrotchAsymUp;
    CrossTypes["CrotchAsymDn"]          = CEPC::kCrotchAsymDn;
    CrossTypes["Legs"]                  = CEPC::kLegs;
    CrossTypes["FlareLegUp"]            = CEPC::kFlareLegUp;
    CrossTypes["FlareLegDn"]           = CEPC::kFlareLegDn;

    std::map < std::string, CEPC::ECrossType>::const_iterator it = CrossTypes.find(type);
    if ( it == CrossTypes.end() ) {
      std::string ms = "Unknown Crossing Type for this geometry " + type;
      throw std::runtime_error(ms.c_str());
    }
    return it->second;
  }
}
#endif // Other_Helpers_hh
