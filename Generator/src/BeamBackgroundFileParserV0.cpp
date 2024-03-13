
#include "BeamBackgroundFileParserV0.h"
#include <sstream>
#include <cmath>

BeamBackgroundFileParserV0::BeamBackgroundFileParserV0(const std::string& filename,
                                                       int pdgid,
                                                       double beam_energy) {
    m_input.open(filename.c_str());
    m_pdgid = pdgid;
    m_beam_energy = beam_energy;
}

bool BeamBackgroundFileParserV0::load(IBeamBackgroundFileParser::BeamBackgroundData& result) {

    if (not m_input.good()) {
        return false;
    }

    // read one record
    std::string tmpline;
    // the format
    double generation_point;
    int loss_turn;
    double z; // unit: m
    double x; // unit: m
    double y; // unit: m
    double cosx; // 
    double cosy; //
    double dz; // unit: m
    double dp; // unit: relative to the E

    while(m_input.good()) {
        std::getline(m_input, tmpline);
        std::stringstream ss;
        ss << tmpline;
        ss >> generation_point; if (ss.fail()) { continue; }
        ss >> loss_turn;        if (ss.fail()) { continue; }
        ss >> z;                if (ss.fail()) { continue; }
        ss >> x;                if (ss.fail()) { continue; }
        ss >> cosx;             if (ss.fail()) { continue; }
        ss >> y;                if (ss.fail()) { continue; }
        ss >> cosy;             if (ss.fail()) { continue; }
        ss >> dz;               if (ss.fail()) { continue; }
        ss >> dp;               if (ss.fail()) { continue; }

        double p = m_beam_energy*(1+dp);

        // Now, we get a almost valid data
        const double m2mm = 1e3; // convert from m to mm
        result.pdgid = m_pdgid;
        result.x     = x * m2mm; 
        result.y     = y * m2mm; 
        result.z     = (z+dz) * m2mm;

        result.px    = p * cosx;
        result.py    = p * cosy;
        result.pz    = p * std::sqrt(1-cosx*cosx-cosy*cosy);

        result.mass  = 0.000511; // assume e-/e+, mass is 0.511 MeV

        return true;
    }
    return false;
}
