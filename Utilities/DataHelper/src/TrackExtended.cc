#include "DataHelper/ClusterExtended.h"
#include "DataHelper/TrackerHitExtended.h"
#include "DataHelper/TrackExtended.h"
#include <math.h>

TrackExtended::TrackExtended( ) {
    _track = NULL;
    _superCluster = NULL;
    _trackerHitVector.clear();
    _clusterVec.clear();
    _group = NULL;
}

TrackExtended::TrackExtended(ConstTrack track) {
    _track = track;
    _superCluster = NULL;
    _trackerHitVector.clear();
    _clusterVec.clear();
    _group = NULL;
    
}

TrackExtended::TrackExtended( TrackerHitExtended * trackerhit) {
    _trackerHitVector.clear();
    _trackerHitVector.push_back(trackerhit);
    _track = NULL;
    _superCluster = NULL;
    _clusterVec.clear();
    _group = NULL;
}


TrackExtended::~TrackExtended() {}

ConstTrack TrackExtended::getTrack() {
  if(!_track.isAvailable()){
    std::cout << "Error: track not available" << _track.isAvailable() << " id= " << _track.id() << std::endl;
    throw std::runtime_error("Error: track not available");
  }
  return _track;
}

const float * TrackExtended::getSeedPosition() {
    return _seedPosition;
}

const float * TrackExtended::getSeedDirection() {
    return _seedDirection;
}

ClusterExtendedVec & TrackExtended::getClusterVec() {
    return _clusterVec;
}

ClusterExtended * TrackExtended::getSuperCluster() {
    return _superCluster;
}

TrackerHitExtendedVec & TrackExtended::getTrackerHitExtendedVec() {
    return _trackerHitVector;
}

void TrackExtended::addCluster(ClusterExtended * cluster) {
    _clusterVec.push_back(cluster);
}

void TrackExtended::setSuperCluster(ClusterExtended * superCluster) {
    _superCluster = superCluster;
}

void TrackExtended::setSeedDirection( float * seedDirection ) {
    _seedDirection[0] = seedDirection[0];
    _seedDirection[1] = seedDirection[1];
    _seedDirection[2] = seedDirection[2];
}

void TrackExtended::setSeedPosition( float * seedPosition) {
    _seedPosition[0] = seedPosition[0];
    _seedPosition[1] = seedPosition[1];
    _seedPosition[2] = seedPosition[2];
}

void TrackExtended::addTrackerHitExtended( TrackerHitExtended * trackerhit) {
  _trackerHitVector.push_back(trackerhit);

}

void TrackExtended::ClearTrackerHitExtendedVec() {
  _trackerHitVector.clear();
}

void TrackExtended::setX0( float x0 ) {
  _x0 = x0;
}

void TrackExtended::setD0( float d0 ) {
  _d0 = d0;
}

void TrackExtended::setZ0( float z0 ) {
  _z0 = z0;
}

void TrackExtended::setY0( float y0 ) {
  _y0 = y0;
}

void TrackExtended::setR0( float r0 ) {
  _r0 = r0;
}

void TrackExtended::setBz( float bz ) {
  _bz = bz;
}

void TrackExtended::setPhi0( float phi0 ) {
  _phi0 = phi0;
}

void TrackExtended::setPhi( float phi ) {
  _phi = phi;
}

void TrackExtended::setOmega( float omega ) {
  _omega = omega;
}

void TrackExtended::setTanLambda( float tanLambda ) {
  _tanLambda = tanLambda;
}
void TrackExtended::setNDF( int ndf) {
  _ndf = ndf;
}


float TrackExtended::getX0() {
  return _x0;
}

float TrackExtended::getY0() {
  return _y0;
}

float TrackExtended::getR0() {
  return _r0;
}

float TrackExtended::getBz() {
  return _bz;
}

float TrackExtended::getD0() {
  return _d0;
}

float TrackExtended::getZ0() {
  return _z0;
}

float TrackExtended::getPhi0() {
  return _phi0;
}

float TrackExtended::getPhi() {
  return _phi;
}

float TrackExtended::getOmega() {
  return _omega;
}

float TrackExtended::getTanLambda() {
  return _tanLambda;
}

int TrackExtended::getNDF() {
  return _ndf;
}

void TrackExtended::setStart(float * xStart) {
  _xStart[0] = xStart[0];
  _xStart[1] = xStart[1];
  _xStart[2] = xStart[2];
}

void TrackExtended::setEnd(float * xEnd) {
  _xEnd[0] = xEnd[0];
  _xEnd[1] = xEnd[1];
  _xEnd[2] = xEnd[2];
}

float * TrackExtended::getStart() {
  return _xStart;
}

float * TrackExtended::getEnd() {
  return _xEnd;
}

void TrackExtended::setGroupTracks( GroupTracks * group ) {
  _group = group;
}

GroupTracks * TrackExtended::getGroupTracks() {
  return _group;
}

float TrackExtended::getChi2() {
  return _chi2;
}

void TrackExtended::setChi2(float chi2) {
  _chi2 = chi2;
}

void TrackExtended::setCovMatrix(float * cov) {
  for (int i=0; i<15; ++i)
    _cov[i] = cov[i];
}

float * TrackExtended::getCovMatrix() {
  return _cov;
}

// void TrackExtended::setType(int type) {
//   _type = type

// }

// int TrackExtended::getType() {
//   return _type;

// }
