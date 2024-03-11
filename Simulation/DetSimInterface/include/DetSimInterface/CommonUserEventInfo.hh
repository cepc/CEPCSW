#ifndef CommonUserEventInfo_hh
#define CommonUserEventInfo_hh

/*
 * Description:
 *   This class is a part of simulation framework to allow users to extend the G4Event.
 *
 *   For example, when G4 converts the EDM4hep/G4 primary vertex/particle to G4 track, 
 *   the relationship between the EDM4hep track and G4 track is missing. 
 *   So a map is used as a bookkeeping.
 *
 * Author:
 *   Tao Lin <lintao AT ihep.ac.cn>
 */

#include "G4VUserEventInformation.hh"
#include <map>

class CommonUserEventInfo: public G4VUserEventInformation {
public:

    CommonUserEventInfo();
    virtual ~CommonUserEventInfo();

public:
    virtual void Print() const;

    // set the relationship between idx in geant4 and idx in mc particle collection.
    // idxG4: G4 track ID (starts from 1)
    // idxEdm4hep: index in MC Particle collection (starts from 0)
    bool setIdxG4Track2Edm4hep(int idxG4, int idxEdm4hep);
    int idxG4Track2Edm4hep(int idxG4) const;
    void dumpIdxG4Track2Edm4hep() const;

private:
    std::map<int, int> m_g4track_to_edm4hep;
};

#endif
