#include "DataHelper/TrackHitPair.h"

TrackHitPair::TrackHitPair(TrackExtended * trkExt, TrackerHitExtended * hitExt, float distance) {
    _trackExtended = trkExt;
    _trackerHitExtended = hitExt;
    _distance = distance;
}
TrackHitPair::~TrackHitPair() {

}

void TrackHitPair::setTrackExtended(TrackExtended * trkExt) {
    _trackExtended = trkExt;
}

void TrackHitPair::setTrackerHitExtended(TrackerHitExtended * hitExt) {
    _trackerHitExtended = hitExt;
}

void TrackHitPair::setDistance(float distance) {
    _distance = distance;
}

TrackExtended * TrackHitPair::getTrackExtended() {
    return _trackExtended;
}

TrackerHitExtended * TrackHitPair::getTrackerHitExtended() {
    return _trackerHitExtended;
}

float TrackHitPair::getDistance() {
    return _distance;
}
