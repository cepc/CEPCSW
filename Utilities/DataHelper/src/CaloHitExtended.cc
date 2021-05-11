#include "DataHelper/CaloHitExtended.h"
#include <math.h>

CaloHitExtended::CaloHitExtended(const edm4hep::CalorimeterHit& calhit, int type)
  :_calohit(calhit){
  _type    = type;
  _index   = 0;
  _calohitTo = NULL;
  _calohitFrom = NULL;
}

CaloHitExtended::~CaloHitExtended(){    
}

CalorimeterHit* CaloHitExtended::getCalorimeterHit() {    
    return &_calohit;
}

CaloHitExtended* CaloHitExtended::getCaloHitTo() {
    return _calohitTo;
}

CaloHitExtended* CaloHitExtended::getCaloHitFrom() {
    return _calohitFrom;
}

ClusterExtended* CaloHitExtended::getClusterExtended() {
    return _clusterAR;    
}

int CaloHitExtended::getIndex() {
    return _index;

}

int CaloHitExtended::getType() {
    return _type;

}
const float* CaloHitExtended::getDirVec() {
    return _dirVec;
}

float  CaloHitExtended::getYresTo() {
    return _yresTo;
}

float  CaloHitExtended::getYresFrom() {
    return _yresFrom;
}



float  CaloHitExtended::getGenericDistance() {
    return _genericDistance;
}

//void CaloHitExtended::setCalorimeterHit(CalorimeterHit* calhit) {
//    _calohit = calhit;
//}

void CaloHitExtended::setCaloHitTo(CaloHitExtended* calhitto) {
    _calohitTo = calhitto; 
}

void CaloHitExtended::setCaloHitFrom(CaloHitExtended* calhitfrom) {
    _calohitFrom = calhitfrom; 
}

void CaloHitExtended::setClusterExtended(ClusterExtended* cluster) {
    _clusterAR = cluster;
}

void CaloHitExtended::setIndex(int index) {
    _index = index;
}

void CaloHitExtended::setType(int type) {
    _type = type;
}

void CaloHitExtended::setDirVec(float* dirVec) {    
    float modulus(0.);
    for (int i(0); i < 3; ++i) {
	modulus += dirVec[i]*dirVec[i];
    }

    modulus = sqrt(modulus);

    if (modulus <= 0.) {
	_dirVec[0] = 0.;
	_dirVec[1] = 0.;
	_dirVec[2] = 1.;
    }
    else {
	_dirVec[0] = dirVec[0] / modulus;
	_dirVec[1] = dirVec[1] / modulus;
	_dirVec[2] = dirVec[2] / modulus;
    }
}

void CaloHitExtended::setYresTo(float yresto) {
    _yresTo = yresto;
}

void CaloHitExtended::setYresFrom(float yresfrom) {
    _yresFrom = yresfrom;
}

void CaloHitExtended::setGenericDistance(float genericDistance) {
    _genericDistance = genericDistance;
}

void CaloHitExtended::setDistanceToCalo(float caloDistance) {
  _caloDistance = caloDistance;
}

float CaloHitExtended::getDistanceToCalo() {
  return _caloDistance;
}

void CaloHitExtended::setDistanceToNearestHit(float distanceToNearest) {
  _distanceToNearestHit = distanceToNearest;
}

float CaloHitExtended::getDistanceToNearestHit() {
  return _distanceToNearestHit;
}
