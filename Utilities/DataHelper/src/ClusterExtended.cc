#include "DataHelper/CaloHitExtended.h"
#include "DataHelper/TrackExtended.h"
#include "DataHelper/ClusterExtended.h"
#include <math.h>

ClusterExtended::ClusterExtended() {
    _hitVector.clear();
    _trackVector.clear();

    _direction[0] = 0.;
    _direction[1] = 0.;
    _direction[2] = 0.;

    _startingPoint[0] = 0.;
    _startingPoint[1] = 0.;
    _startingPoint[2] = 0.;	

}

ClusterExtended::ClusterExtended( Cluster * cluster ) {
    _hitVector.clear();
    _trackVector.clear();

    _direction[0] = 0.;
    _direction[1] = 0.;
    _direction[2] = 0.;

    _startingPoint[0] = 0.;
    _startingPoint[1] = 0.;
    _startingPoint[2] = 0.;	

    _cluster = cluster;

}




ClusterExtended::ClusterExtended(CaloHitExtended * calohit) {

    _hitVector.clear();
    _hitVector.push_back(calohit);
    _trackVector.clear();

    float rad(0);

    _startingPoint[0] = calohit->getCalorimeterHit()->getPosition().x;
    _startingPoint[1] = calohit->getCalorimeterHit()->getPosition().y;
    _startingPoint[2] = calohit->getCalorimeterHit()->getPosition().z;
    
    for (int i(0); i < 3; ++i) {
	rad += _startingPoint[i]*_startingPoint[i];
    }

    rad = sqrt(rad);

    _direction[0] = _startingPoint[0]/rad;
    _direction[1] = _startingPoint[1]/rad;
    _direction[2] = _startingPoint[2]/rad;
    
}

ClusterExtended::ClusterExtended(TrackExtended * track) {
    _hitVector.clear();
    _trackVector.clear();
    _trackVector.push_back(track);
    for (int i(0); i < 3; ++i) {
	_startingPoint[i] = track->getSeedPosition()[i];
	_direction[i] = track->getSeedDirection()[i];
    }
}


ClusterExtended::~ClusterExtended() {
   _hitVector.clear(); 
   _trackVector.clear();
}

CaloHitExtendedVec& ClusterExtended::getCaloHitExtendedVec() {
    return _hitVector;
}

TrackExtendedVec& ClusterExtended::getTrackExtendedVec() {
    return _trackVector;
}

const float* ClusterExtended::getStartingPoint() {
    return _startingPoint;
}

const float* ClusterExtended::getDirection() {
    return _direction;
}

void ClusterExtended::setStartingPoint(float* sPoint) {
    _startingPoint[0] = sPoint[0];
    _startingPoint[1] = sPoint[1];
    _startingPoint[2] = sPoint[2];
    
}

void ClusterExtended::setDirection(float* direct) {
    _direction[0] = direct[0];
    _direction[1] = direct[1];
    _direction[2] = direct[2];
}

void ClusterExtended::addCaloHitExtended(CaloHitExtended* calohit) {
    _hitVector.push_back(calohit);
}

void ClusterExtended::addTrackExtended(TrackExtended * track) {
    _trackVector.push_back(track);
}

void ClusterExtended::Clear() {
    _hitVector.clear();
    _trackVector.clear();

}

void ClusterExtended::setType( int type ) {
  _type = type;
}

int ClusterExtended::getType() {
  return _type;
}

void ClusterExtended::setCluster(Cluster * cluster) {
  _cluster = cluster;
}

Cluster * ClusterExtended::getCluster() {
  return _cluster;
}

void ClusterExtended::setAxis(float * axis) {
  _axis[0] = axis[0];
  _axis[1] = axis[1];
  _axis[2] = axis[2];
}
float * ClusterExtended::getAxis() {
  return _axis;
}

void ClusterExtended::setEccentricity( float eccentricity) {
  _eccentricity = eccentricity;
}
float ClusterExtended::getEccentricity() {
  return _eccentricity;
}

void ClusterExtended::setHelix(HelixClass helix) {
  _helix = helix;
  int nHits = int(_hitVector.size());
  float timeMax = -1.0e+20;
  float timeMin = 1.0e+20;
  for (int ihit=0;ihit<nHits;++ihit) {
    float pos[3];
    for (int i=0;i<3;++i)
      pos[i]=_hitVector[ihit]->getCalorimeterHit()->getPosition()[i];
      
    float distance[3];
    float time = _helix.getDistanceToPoint(pos,distance);
    if (time > timeMax) {
      timeMax = time;
      _upEdge[0] = pos[0];
      _upEdge[1] = pos[1];
      _upEdge[2] = pos[2];
    }
    if (time < timeMin) {
      timeMin = time;
      _lowEdge[0] = pos[0];
      _lowEdge[1] = pos[1];
      _lowEdge[2] = pos[2];
    }
      
  }

}
HelixClass & ClusterExtended::getHelix() {
  return _helix;
}

void ClusterExtended::setHelixChi2R(float helixChi2) {
  _helixChi2R = helixChi2;
}
float ClusterExtended::getHelixChi2R() {
  return _helixChi2R;
}

void ClusterExtended::setHelixChi2Z(float helixChi2) {
  _helixChi2Z = helixChi2;
}
float ClusterExtended::getHelixChi2Z() {
  return _helixChi2Z;
}



void ClusterExtended::setPosition(float * position) {
  _position[0] = position[0];
  _position[1] = position[1];
  _position[2] = position[2];

}

float * ClusterExtended::getPosition() {
  return _position;
}

void ClusterExtended::setUpEdge(float * upEdge) {
  _upEdge[0] = upEdge[0];
  _upEdge[1] = upEdge[1];
  _upEdge[2] = upEdge[2];
}

void ClusterExtended::setLowEdge(float * lowEdge) {
  _lowEdge[0] = lowEdge[0]; 
  _lowEdge[1] = lowEdge[1]; 
  _lowEdge[2] = lowEdge[2]; 
}

float * ClusterExtended::getUpEdge() {
  return _upEdge;
}

float * ClusterExtended::getLowEdge() {
  return _lowEdge;
}

