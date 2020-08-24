#include <ILDImpl/MiniVector.h>

#include "TrackSystemSvc/IMarlinTrack.h"
#include "CLHEP/Vector/ThreeVector.h"

using namespace KiTrack;
using namespace KiTrackMarlin;

MiniVector::MiniVector(edm4hep::TrackerHit * outer, edm4hep::TrackerHit * inner) { 
  HitVec.push_back(outer);
  HitVec.push_back(inner);
}


MiniVector::MiniVector( TrackerHitVec hitPair){

  HitVec = hitPair ;

}

MiniVector::~MiniVector(){

}


double * MiniVector::getPosition(){

  double x_outer = HitVec[0]->getPosition()[0];
  double y_outer = HitVec[0]->getPosition()[1];
  double z_outer = HitVec[0]->getPosition()[2];
  
  double x_inner = HitVec[1]->getPosition()[0];
  double y_inner = HitVec[1]->getPosition()[1];
  double z_inner = HitVec[1]->getPosition()[2];

  double *_pos = new double[3] ;

  //_pos[0] = (x_outer - x_inner)/2. ;
  //_pos[1] = (y_outer - y_inner)/2. ;
  //_pos[2] = (z_outer - z_inner)/2. ;

  _pos[0] = x_inner ;
  _pos[1] = y_inner ;
  _pos[2] = z_inner ;

  return _pos ;

}


TrackerHitVec MiniVector::getTrackerHitVec(){

  return HitVec ;

}



double MiniVector::getPhi(){

  double x_outer = HitVec[0]->getPosition()[0];
  double y_outer = HitVec[0]->getPosition()[1];
  
  double x_inner = HitVec[1]->getPosition()[0];
  double y_inner = HitVec[1]->getPosition()[1];

  double phiMV =  atan2((y_outer-y_inner),(x_outer-x_inner))  ;

  //std::cout << " outer x,y " << x_outer << " , " << y_outer << " inner x,y " << x_inner << " , " << y_inner << std::endl ;

  //std::cout << " Calling MiniVector::getPhi, returning " << phiMV << std::endl ;

  return phiMV ;

}

double MiniVector::getTheta(){

  double x_outer = HitVec[0]->getPosition()[0];
  double y_outer = HitVec[0]->getPosition()[1];
  double z_outer = HitVec[0]->getPosition()[2];
  
  double x_inner = HitVec[1]->getPosition()[0];
  double y_inner = HitVec[1]->getPosition()[1];
  double z_inner = HitVec[1]->getPosition()[2];

  double thetaMV = atan2(sqrt((x_outer-x_inner)*(x_outer-x_inner) + (y_outer-y_inner)*(y_outer-y_inner)),(z_outer-z_inner)) ;

  return thetaMV ;

}


double * MiniVector::getXYZ(){

  double x_outer = HitVec[0]->getPosition()[0];
  double y_outer = HitVec[0]->getPosition()[1];
  double z_outer = HitVec[0]->getPosition()[2];
  
  double x_inner = HitVec[1]->getPosition()[0];
  double y_inner = HitVec[1]->getPosition()[1];
  double z_inner = HitVec[1]->getPosition()[2];

  double *xyz = new double[3] ;

  xyz[0] = x_outer - x_inner ;
  xyz[1] = y_outer - y_inner ;
  xyz[2] = z_outer - z_inner ;


  return xyz ;

}


double MiniVector::get3DAngleMV(MiniVector *MinVec2){

  double *xyz1 = new double[3] ;
  xyz1 = this->getXYZ() ;
  CLHEP::Hep3Vector v1( xyz1[0], xyz1[1], xyz1[2] ) ;

  double *xyz2 = new double[3] ;
  xyz2 = MinVec2->getXYZ() ;

  CLHEP::Hep3Vector v2( xyz2[0], xyz2[1], xyz2[2] ) ;

  double ScalarProd = v1.dot(v2);
  double magV1 = v1.r();
  double magV2 = v2.r();

  std::cout << "####### MINIVECTOR::get3DAngle, Scalar product " << ScalarProd << " mag. of the instantiate mv " << magV1 << " mag. of the second mv " << magV2 << std::endl ; 

  double Angle3D =  acos(ScalarProd / (magV1*magV2)) ; 
 
  return Angle3D ;

}


