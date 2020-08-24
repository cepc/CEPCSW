#include "DataHelper/TrackExtended.h"
#include "DataHelper/GroupTracks.h"

GroupTracks::GroupTracks() {
  _trackARVec.clear();
  _edges[0] = 0.0;
  _edges[1] = 0.0;
}

GroupTracks::GroupTracks( TrackExtended * track ) {
  _trackARVec.clear();
  _trackARVec.push_back( track );
  _edges[0] = 0.0;
  _edges[1] = 0.0;
}

GroupTracks::~GroupTracks() {}

void GroupTracks::addTrackExtended( TrackExtended * track ) {  
  _trackARVec.push_back( track );
}

void GroupTracks::ClearTrackExtendedVec() {
  _trackARVec.clear();
}

TrackExtendedVec & GroupTracks::getTrackExtendedVec() {
  return _trackARVec;
}

void GroupTracks::setEdges(float * edges) {

  _edges[0] = edges[0];
  _edges[1] = edges[1];  

}

float * GroupTracks::getEdges() {
  return _edges;
}
