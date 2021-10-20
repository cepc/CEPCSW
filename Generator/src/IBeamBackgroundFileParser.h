#ifndef IBeamBackgroundFileParser_h
#define IBeamBackgroundFileParser_h

/*
 * Description:
 *   This interface is used to load the beam background information, such as:
 *     - pdgid (optional)
 *         About the pdgid, it will be e+/e- in most cases.
 *     - x/y/z
 *     - t (optional) 
 *     - px/py/pz
 *         About the time, it could be set in the GtBeamBackgroundTool.
 *
 * Author:
 *   Tao Lin <lintao AT ihep.ac.cn>
 */

class IBeamBackgroundFileParser {
public:
    // Internal used Data
    struct BeamBackgroundData {
        int pdgid;

        double x; // unit: mm
        double y; // unit: mm
        double z; // unit: mm
        double t; // unit: ns

        double px; // unit: GeV
        double py; // unit: GeV
        double pz; // unit: GeV
        double mass; // unit: GeV

        BeamBackgroundData() 
          : pdgid(11), x(0), y(0), z(0), t(0), 
            px(0), py(0), pz(0), mass(0) {}
        
    };

    // return false if failed to load the data
    virtual bool load(BeamBackgroundData&) = 0;
};

#endif
