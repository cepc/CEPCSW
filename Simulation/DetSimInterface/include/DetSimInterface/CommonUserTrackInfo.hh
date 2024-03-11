#ifndef CommonUserTrackInfo_hh
#define CommonUserTrackInfo_hh

/* Description:
 *   This class is a part of simulation framework to extend the G4Track.
 *
 *   Some secondaries are created due to decay. However, their G4 Track IDs are
 *   not valid until the tracks are tracking by geant4. In order to associate
 *   these tracks and their edm4hep MC particle, we use the track information
 *   to record the extra track information.
 *
 * Author:
 *   Tao Lin <lintao AT ihep.ac.cn>
 */

#include "G4VUserTrackInformation.hh"

class CommonUserTrackInfo: public G4VUserTrackInformation {
public:
    CommonUserTrackInfo();
    ~CommonUserTrackInfo();

public:

    virtual void Print() const;

    // get the idx in the EDM4hep MC particle collection
    bool setIdxEdm4hep(int idxEdm4hep);
    int idxEdm4hep() const;

private:
    int m_idxEdm4hep = -1;
};

#endif
