#include "DataHelper/TrackExtended.h"
#include "DataHelper/TrackerHitExtended.h"

TrackerHitExtended::TrackerHitExtended(const edm4hep::ConstTrackerHit trackerhit):
  _trackerHit(trackerhit){
  _trackAR = NULL;
  _trackVecAR.clear();
  _usedInFit = false;
}

TrackerHitExtended::~TrackerHitExtended() {

}

void TrackerHitExtended::setTrackExtended(TrackExtended * trackAR) {
    _trackAR = trackAR;
    _trackVecAR.clear();
    _trackVecAR.push_back(trackAR);
      
} 

void TrackerHitExtended::addTrackExtended(TrackExtended * trackAR) {
  _trackVecAR.push_back(trackAR);
}

TrackExtendedVec & TrackerHitExtended::getTrackExtendedVec() {
  return _trackVecAR;
}


void TrackerHitExtended::setTrackerHitTo(TrackerHitExtended * hitTo) {
    _hitTo = hitTo;
}

void TrackerHitExtended::setTrackerHitFrom(TrackerHitExtended * hitFrom) {
    _hitFrom = hitFrom;

} 

void TrackerHitExtended::setGenericDistance(float genericDistance) {
    _genericDistance = genericDistance; 
}

//void TrackerHitExtended::setTrackerHit(edm4hep::ConstTrackerHit hit) {
//    _trackerHit = hit;
//}

void TrackerHitExtended::setYresTo(float yresTo) {
    _yresTo = yresTo;
}

void TrackerHitExtended::setYresFrom(float yresFrom) {
    _yresFrom = yresFrom;
}

void TrackerHitExtended::setDirVec(float * dirVec) {
    _dirVec[0] = dirVec[0];
    _dirVec[1] = dirVec[1];
    _dirVec[2] = dirVec[2];
}

void TrackerHitExtended::setType(int ityp) {
    _type = ityp;
}

void TrackerHitExtended::setDet(int idet) {
    _idet = idet;
}

edm4hep::ConstTrackerHit TrackerHitExtended::getTrackerHit() {
    return _trackerHit;
}

TrackExtended * TrackerHitExtended::getTrackExtended() {
    return _trackAR;
}

TrackerHitExtended * TrackerHitExtended::getTrackerHitFrom() {
    return _hitFrom;
}
TrackerHitExtended * TrackerHitExtended::getTrackerHitTo() {
    return _hitTo;
}
float TrackerHitExtended::getGenericDistance() {
    return _genericDistance;
}
float TrackerHitExtended::getYresTo() {
    return _yresTo;
}
float TrackerHitExtended::getYresFrom() {
    return _yresFrom;
}

float * TrackerHitExtended::getDirVec() {
    return _dirVec;
}

void TrackerHitExtended::clearTrackVec() {
  _trackVecAR.clear();
}

float TrackerHitExtended::getResolutionRPhi() {
  return _rphiReso;
}

float TrackerHitExtended::getResolutionZ() {
  return _zReso;
}

void TrackerHitExtended::setResolutionRPhi(float rphiReso) {
  _rphiReso = rphiReso;
}

void TrackerHitExtended::setResolutionZ(float zReso) {
  _zReso = zReso;
}

int TrackerHitExtended::getType() {
    return _type;
}

int TrackerHitExtended::getDet() {
    return _idet;
}

void TrackerHitExtended::setUsedInFit(bool usedInFit) {
    _usedInFit = usedInFit;
}

bool TrackerHitExtended::getUsedInFit() {
    return _usedInFit;
}

