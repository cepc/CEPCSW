#include "Utility.h"

std::string Convert (float number){
    std::ostringstream buff;
    buff<<number;
    return buff.str();   
}
