#ifndef BeamBackgroundFileParserV0_h
#define BeamBackgroundFileParserV0_h

#include "IBeamBackgroundFileParser.h"

#include <fstream>

class BeamBackgroundFileParserV0: public IBeamBackgroundFileParser {
public:
    BeamBackgroundFileParserV0(const std::string& filename);

    bool load(IBeamBackgroundFileParser::BeamBackgroundData&);
private:
    std::ifstream m_input;
};

#endif
