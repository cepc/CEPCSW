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

class CommonUserEventInfo: public G4VUserEventInformation {
public:

    CommonUserEventInfo();
    virtual ~CommonUserEventInfo();

public:
    virtual void Print() const;

};

#endif
