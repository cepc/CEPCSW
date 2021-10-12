
#include "BeamBackgroundFileParserV0.h"

BeamBackgroundFileParserV0::BeamBackgroundFileParserV0(const std::string& filename) {
    m_input.open(filename.c_str());
}

bool BeamBackgroundFileParserV0::load(IBeamBackgroundFileParser::BeamBackgroundData& result) {
    return true;
}
