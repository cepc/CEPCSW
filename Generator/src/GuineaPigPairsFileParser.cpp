#include "GuineaPigPairsFileParser.h"
#include <sstream>
#include <cmath>

GuineaPigPairsFileParser::GuineaPigPairsFileParser(const std::string& filename) {
    m_input.open(filename.c_str());
}

bool GuineaPigPairsFileParser::load(IBeamBackgroundFileParser::BeamBackgroundData& result) {
    if (not m_input.good()) {
        return false;
    }

    // read one record
    std::string tmpline;
    // the format
    double energy; // unit: GeV
    double vx; // unit: c
    double vy; // unit: c
    double vz; // unit: c
    double x; // unit: nm
    double y; // unit: nm
    double z; // unit: nm
    int process;

    while(m_input.good()) {
        std::getline(m_input, tmpline);
        std::stringstream ss;
        ss << tmpline;
        ss >> energy;           if (ss.fail()) { continue; }
        ss >> vx;               if (ss.fail()) { continue; }
        ss >> vy;               if (ss.fail()) { continue; }
        ss >> vz;               if (ss.fail()) { continue; }
        ss >> x;                if (ss.fail()) { continue; }
        ss >> y;                if (ss.fail()) { continue; }
        ss >> z;                if (ss.fail()) { continue; }
        ss >> process;          if (ss.fail()) { continue; }

        int pdgid = 11; // 11: electron; -11: positron
        if (energy<0) pdgid = -11;

        double p = std::fabs(energy);
        double v = sqrt(vx*vx+vy*vy+vz*vz);

        // Now, we get a almost valid data
        const double nm2mm = 1e-6; // convert from nm to mm
        result.pdgid = pdgid;
        result.x     = x * nm2mm; 
        result.y     = y * nm2mm; 
        result.z     = z * nm2mm;

        result.px    = p * vx/v;
        result.py    = p * vy/v;
        result.pz    = p * vz/v;

        result.mass  = 0.000511; // assume e-/e+, mass is 0.511 MeV

        return true;
    }
    return false;

}
