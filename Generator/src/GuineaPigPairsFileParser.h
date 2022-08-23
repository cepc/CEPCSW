#ifndef GuineaPigPairsFileParser_h
#define GuineaPigPairsFileParser_h

#include "IBeamBackgroundFileParser.h"
#include <fstream>

/* Format of Guinea-Pig Pairs:
 *
 *   E vx vy vz x y z process
 *
 * Notes:
 *   - E (GeV). If E>0, it is electron. If E<0, it is positron
 *   - vx/vy/vz (speed of light)
 *   - x/y/z (nm)
 *   - process
 *     - 0: Breit-Wheeler
 *     - 1: Bethe-Heitler
 *     - 2: Landau-Lifschitz
 *
 */

class GuineaPigPairsFileParser: public IBeamBackgroundFileParser {
public:
    GuineaPigPairsFileParser(const std::string& filename);

    bool load(IBeamBackgroundFileParser::BeamBackgroundData&);

private:
    std::ifstream m_input;
};

#endif
