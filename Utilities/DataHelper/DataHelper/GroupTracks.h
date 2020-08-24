#ifndef GROUPTRACKS_H
#define GROUPTRACKS_H 1

//#include "TrackExtended.h"
//#include "ClusterExtended.h"
#include <vector>

//fg : forwar declaration needed because of circular include ....
class TrackExtended;
typedef std::vector<TrackExtended*> TrackExtendedVec;
//fg : forward ....

/**
 * Class GroupTracks is needed to group track segments <br>
 * with close helix parameters. Needed for Tracking. <br>
 *  * @author A. Raspereza (DESY)<br>
 * @version $Id: GroupTracks.h,v 1.4 2007-09-05 09:39:49 rasp Exp $<br>
 */
//class GroupTracks;

//typedef std::vector<GroupTracks*>  GroupTracksVec;

class GroupTracks {

 public:
  GroupTracks();
  GroupTracks(TrackExtended * track );
  ~GroupTracks();

  void addTrackExtended( TrackExtended * track );
  void ClearTrackExtendedVec();
  TrackExtendedVec & getTrackExtendedVec();
  void setEdges(float * edges);
  float * getEdges();

 private:

  TrackExtendedVec _trackARVec;

  float _edges[2];

};

typedef std::vector<GroupTracks*>  GroupTracksVec;

#endif
