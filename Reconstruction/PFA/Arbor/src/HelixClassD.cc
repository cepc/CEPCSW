#include "HelixClassD.hh"
#include <math.h>
#include <stdlib.h>
#include <iostream>
//#include "ced_cli.h"

HelixClassD::HelixClassD() {
    _const_2pi = 2.0*M_PI;
    _const_pi2 = 0.5*M_PI;
    _FCT = 2.99792458E-4;
}

HelixClassD::~HelixClassD() {}

void HelixClassD::Initialize_VP(float * pos, float * mom, float q, float B) {
    _referencePoint[0] = pos[0];
    _referencePoint[1] = pos[1];
    _referencePoint[2] = pos[2];
    _momentum[0] = mom[0];
    _momentum[1] = mom[1];
    _momentum[2] = mom[2];
    _charge = q;
    _bField = B;
    _pxy = sqrt(mom[0]*mom[0]+mom[1]*mom[1]);
    _radius = _pxy / (_FCT*B);
    _omega = q/_radius;
    _tanLambda = mom[2]/_pxy;
    _phiMomRefPoint = atan2(mom[1],mom[0]);
    _xCentre = pos[0] + _radius*cos(_phiMomRefPoint-_const_pi2*q);
    _yCentre = pos[1] + _radius*sin(_phiMomRefPoint-_const_pi2*q);
    _phiRefPoint = atan2(pos[1]-_yCentre,pos[0]-_xCentre);
    _phiAtPCA = atan2(-_yCentre,-_xCentre);
    _phi0 = -_const_pi2*q + _phiAtPCA;
    while (_phi0<0) _phi0+=_const_2pi;
    while (_phi0>=_const_2pi) _phi0-=_const_2pi;
    _xAtPCA = _xCentre + _radius*cos(_phiAtPCA);
    _yAtPCA = _yCentre + _radius*sin(_phiAtPCA);
    //    _d0 = -_xAtPCA*sin(_phi0) + _yAtPCA*cos(_phi0);
    double pxy = double(_pxy);
    double radius = pxy/double(_FCT*B);
    double xCentre = double(pos[0]) + radius*double(cos(_phiMomRefPoint-_const_pi2*q));
    double yCentre = double(pos[1]) + radius*double(sin(_phiMomRefPoint-_const_pi2*q));
    
    double d0;

    if (q>0) {
      d0 = double(q)*radius - double(sqrt(xCentre*xCentre+yCentre*yCentre));
    }
    else {
      d0 = double(q)*radius + double(sqrt(xCentre*xCentre+yCentre*yCentre));
    }

    _d0 = float(d0);

//     if (fabs(_d0)>0.001 ) {
//       std::cout << "New helix : " << std::endl;
//       std::cout << " Position : " << pos[0] 
// 		<< " " << pos[1]
// 		<< " " << pos[2] << std::endl;
//       std::cout << " Radius = " << _radius << std::endl;
//       std::cout << " RC = " << sqrt(_xCentre*_xCentre+_yCentre*_yCentre) << std::endl;  
//       std::cout << " D0 = " << _d0 << std::endl;
//     }

    _pxAtPCA = _pxy*cos(_phi0);
    _pyAtPCA = _pxy*sin(_phi0);
    float deltaPhi = _phiRefPoint - _phiAtPCA;    
    float xCircles = -pos[2]*q/(_radius*_tanLambda) - deltaPhi;
    xCircles = xCircles/_const_2pi;
    int nCircles;
    int n1,n2;

    if (xCircles >= 0.) {
	n1 = int(xCircles);
	n2 = n1 + 1;
    }
    else {
	n1 = int(xCircles) - 1;
	n2 = n1 + 1;
    }
    
    if (fabs(n1-xCircles) < fabs(n2-xCircles)) {
	nCircles = n1;
    }
    else {
	nCircles = n2;
    }
    _z0 = pos[2] + _radius*_tanLambda*q*(deltaPhi + _const_2pi*nCircles);

}

void HelixClassD::Initialize_Canonical(float phi0, float d0, float z0, 
				      float omega, float tanLambda, float B) {
    _omega = omega;
    _d0 = d0;
    _phi0 = phi0;
    _z0 = z0;
    _tanLambda = tanLambda;
    _charge = omega/fabs(omega);
    _radius = 1./fabs(omega);
    _xAtPCA = -_d0*sin(_phi0);
    _yAtPCA = _d0*cos(_phi0);    
    _referencePoint[0] = _xAtPCA;
    _referencePoint[1] = _yAtPCA;
    _referencePoint[2] = _z0;
    _pxy = _FCT*B*_radius;
    _momentum[0] = _pxy*cos(_phi0);
    _momentum[1] = _pxy*sin(_phi0);
    _momentum[2] = _tanLambda * _pxy;    
    _pxAtPCA = _momentum[0];
    _pyAtPCA = _momentum[1];
    _phiMomRefPoint = atan2(_momentum[1],_momentum[0]);
    _xCentre = _referencePoint[0] + 
      _radius*cos(_phi0-_const_pi2*_charge);
    _yCentre = _referencePoint[1] + 
      _radius*sin(_phi0-_const_pi2*_charge);
    _phiAtPCA = atan2(-_yCentre,-_xCentre);
    _phiRefPoint =  _phiAtPCA ;
    _bField = B;
}


void HelixClassD::Initialize_BZ(float xCentre, float yCentre, float radius, 
			       float bZ, float phi0, float B, float signPz,
			       float zBegin) {

  _radius = radius;
  _pxy = _FCT*B*_radius;
  _charge = -(bZ*signPz)/fabs(bZ*signPz);
  _momentum[2] = -_charge*_pxy/(bZ*_radius);
  _xCentre = xCentre;
  _yCentre = yCentre;
  _omega = _charge/radius;
  _phiAtPCA = atan2(-_yCentre,-_xCentre);
  _phi0 = -_const_pi2*_charge + _phiAtPCA;
  while (_phi0<0) _phi0+=_const_2pi;
  while (_phi0>=_const_2pi) _phi0-=_const_2pi;
  _xAtPCA = _xCentre + _radius*cos(_phiAtPCA);
  _yAtPCA = _yCentre + _radius*sin(_phiAtPCA);
  _d0 = -_xAtPCA*sin(_phi0) + _yAtPCA*cos(_phi0);
  _pxAtPCA = _pxy*cos(_phi0);
  _pyAtPCA = _pxy*sin(_phi0);
  _referencePoint[2] = zBegin;
  _referencePoint[0] = xCentre + radius*cos(bZ*zBegin+phi0);
  _referencePoint[1] = yCentre + radius*sin(bZ*zBegin+phi0);
  _phiRefPoint = atan2(_referencePoint[1]-_yCentre,_referencePoint[0]-_xCentre);
  _phiMomRefPoint =  -_const_pi2*_charge + _phiRefPoint;
  _tanLambda = _momentum[2]/_pxy;
  _momentum[0] = _pxy*cos(_phiMomRefPoint);
  _momentum[1] = _pxy*sin(_phiMomRefPoint);
  
  float deltaPhi = _phiRefPoint - _phiAtPCA;    
  float xCircles = bZ*_referencePoint[2] - deltaPhi;
  xCircles = xCircles/_const_2pi;
  int nCircles;
  int n1,n2;

  if (xCircles >= 0.) {
    n1 = int(xCircles);
    n2 = n1 + 1;
  }
  else {
    n1 = int(xCircles) - 1;
    n2 = n1 + 1;
  }
  
  if (fabs(n1-xCircles) < fabs(n2-xCircles)) {
    nCircles = n1;
  }
  else {
    nCircles = n2;
  }  
  _z0 = _referencePoint[2] - (deltaPhi + _const_2pi*nCircles)/bZ;  
  _bField = B;

}

const float * HelixClassD::getMomentum() {
    return _momentum;
}
const float * HelixClassD::getReferencePoint() {
    return _referencePoint;
}
float HelixClassD::getPhi0() {
  if (_phi0<0.0)
    _phi0 += 2*M_PI;
  return _phi0;
}
float HelixClassD::getD0() {
    return _d0;
}
float HelixClassD::getZ0() {
    return _z0;
}
float HelixClassD::getOmega() {
    return _omega;
}
float HelixClassD::getTanLambda() {
    return _tanLambda;
}
float HelixClassD::getPXY() {
    return _pxy;
}
float HelixClassD::getXC() {
  return _xCentre;
}

float HelixClassD::getYC() {
  return _yCentre;
}

float HelixClassD::getRadius() {
  return _radius;
}

float HelixClassD::getBz() {
  return _bZ;
}

float HelixClassD::getPhiZ() {
  return _phiZ;
}

float HelixClassD::getCharge() {
    return _charge;
}

float HelixClassD::getPointInXY(float x0, float y0, float ax, float ay, 
			      float * ref , float * point) {

  float time;

  float AA = sqrt(ax*ax+ay*ay);


  if (AA <= 0) {
    time = -1.0e+20; 
    return time;
  }


  float BB = ax*(x0-_xCentre) + ay*(y0-_yCentre);
  BB = BB / AA;

  float CC = (x0-_xCentre)*(x0-_xCentre) 
    + (y0-_yCentre)*(y0-_yCentre) - _radius*_radius;

  CC = CC / AA;

  float DET = BB*BB - CC;
  float tt1 = 0.;
  float tt2 = 0.;
  float xx1,xx2,yy1,yy2; 


  if (DET < 0 ) {
    time = 1.0e+10;
    point[0]=0.0;
    point[1]=0.0;
    point[2]=0.0;
    return time;
  }
  
  
  tt1 = - BB + sqrt(DET);
  tt2 = - BB - sqrt(DET);

  xx1 = x0 + tt1*ax;
  yy1 = y0 + tt1*ay;
  xx2 = x0 + tt2*ax;
  yy2 = y0 + tt2*ay;
  
  float phi1 = atan2(yy1-_yCentre,xx1-_xCentre);
  float phi2 = atan2(yy2-_yCentre,xx2-_xCentre);
  float phi0 = atan2(ref[1]-_yCentre,ref[0]-_xCentre);

  float dphi1 = phi1 - phi0;
  float dphi2 = phi2 - phi0;

  if (dphi1 < 0 && _charge < 0) {
    dphi1 = dphi1 + _const_2pi;
  }
  else if (dphi1 > 0 && _charge > 0) { 
    dphi1 = dphi1 - _const_2pi;
  }

  if (dphi2 < 0 && _charge < 0) {
    dphi2 = dphi2 + _const_2pi;
  }
  else if (dphi2 > 0 && _charge > 0) { 
    dphi2 = dphi2 - _const_2pi;
  }

  // Times
  tt1 = -_charge*dphi1*_radius/_pxy;
  tt2 = -_charge*dphi2*_radius/_pxy;

  if (tt1 < 0. )
    std::cout << "WARNING " << tt1 << std::endl;
  if (tt2 < 0. )
    std::cout << "WARNING " << tt2 << std::endl;
  

  if (tt1 < tt2) {
    point[0] = xx1;
    point[1] = yy1;
    time = tt1;
  }
  else {
    point[0] = xx2;
    point[1] = yy2;
    time = tt2;
  }

  point[2] = ref[2]+time*_momentum[2];

  

  return time;

}


float HelixClassD::getPointOnCircle(float Radius, float * ref, float * point) {

  float distCenterToIP = sqrt(_xCentre*_xCentre + _yCentre*_yCentre);

  point[0] = 0.0;
  point[1] = 0.0;
  point[2] = 0.0;

  if ((distCenterToIP+_radius)<Radius) {
    float xx = -1.0e+20;
    return xx;
  }

  if ((_radius+Radius)<distCenterToIP) {
    float xx = -1.0e+20;
    return xx;
  }

  float phiCentre = atan2(_yCentre,_xCentre);
  float phiStar   = Radius*Radius + distCenterToIP*distCenterToIP 
                                    - _radius*_radius;

  phiStar = 0.5*phiStar/fmax(1.0e-20,Radius*distCenterToIP);
  
  if (phiStar > 1.0) 
    phiStar = 0.9999999;
  
  if (phiStar < -1.0)
    phiStar = -0.9999999;
  
  phiStar = acos(phiStar);

  float tt1,tt2,time;

  float xx1 = Radius*cos(phiCentre+phiStar);
  float yy1 = Radius*sin(phiCentre+phiStar);

  float xx2 = Radius*cos(phiCentre-phiStar);
  float yy2 = Radius*sin(phiCentre-phiStar);


  float phi1 = atan2(yy1-_yCentre,xx1-_xCentre);
  float phi2 = atan2(yy2-_yCentre,xx2-_xCentre);
  float phi0 = atan2(ref[1]-_yCentre,ref[0]-_xCentre);

  float dphi1 = phi1 - phi0;
  float dphi2 = phi2 - phi0;

  if (dphi1 < 0 && _charge < 0) {
    dphi1 = dphi1 + _const_2pi;
  }
  else if (dphi1 > 0 && _charge > 0) { 
    dphi1 = dphi1 - _const_2pi;
  }

  if (dphi2 < 0 && _charge < 0) {
    dphi2 = dphi2 + _const_2pi;
  }
  else if (dphi2 > 0 && _charge > 0) { 
    dphi2 = dphi2 - _const_2pi;
  }

  // Times
  tt1 = -_charge*dphi1*_radius/_pxy;
  tt2 = -_charge*dphi2*_radius/_pxy;

  if (tt1 < 0. )
    std::cout << "WARNING " << tt1 << std::endl;
  if (tt2 < 0. )
    std::cout << "WARNING " << tt2 << std::endl;
  

  float time2;
  if (tt1 < tt2) {
    point[0] = xx1;
    point[1] = yy1;
    point[3] = xx2;
    point[4] = yy2;
    time = tt1;
    time2 = tt2;
  }
  else {
    point[0] = xx2;
    point[1] = yy2;
    point[3] = xx1;
    point[4] = yy1;
    time = tt2;
    time2 = tt1;
  }

  point[2] = ref[2]+time*_momentum[2];
  point[5] = ref[2]+time2*_momentum[2];
  

  return time;

}


float HelixClassD::getPointInZ(float zLine, float * ref, float * point) {

  float time = zLine - ref[2];

  if (_momentum[2] == 0.) {
    time = -1.0e+20;
    point[0] = 0.;
    point[1] = 0.;
    point[2] = 0.;
    return time;
  }

  time = time/_momentum[2];

  float phi0 = atan2(ref[1] - _yCentre , ref[0] - _xCentre);
  float phi = phi0 - _charge*_pxy*time/_radius;
  float xx = _xCentre + _radius*cos(phi);
  float yy = _yCentre + _radius*sin(phi);

  point[0] = xx;
  point[1] = yy;
  point[2] = zLine;

  return time;


}

int HelixClassD::getNCircle(float * xPoint) {

  int nCircles = 0;

  float phi = atan2(xPoint[1]-_yCentre,xPoint[0]-_xCentre);
  float phi0 = atan2(_referencePoint[1]-_yCentre,_referencePoint[0]-_xCentre);

  if (fabs(_tanLambda*_radius)>1.0e-20) {
    float xCircles = phi0 - phi -_charge*(xPoint[2]-_referencePoint[2])/(_tanLambda*_radius);
    xCircles = xCircles/_const_2pi;
    int n1,n2;

    if (xCircles >= 0.) {
        n1 = int(xCircles);
        n2 = n1 + 1;
    }
    else {
        n1 = int(xCircles) - 1;
        n2 = n1 + 1;
    }

    if (fabs(n1-xCircles) < fabs(n2-xCircles)) {
        nCircles = n1;
    }
    else {
        nCircles = n2;
    }

  }
  return nCircles;

}

float HelixClassD::getDistanceToPoint(float * xPoint, float * Distance) {

  float zOnHelix;
  float phi = atan2(xPoint[1]-_yCentre,xPoint[0]-_xCentre);
  float phi0 = atan2(_referencePoint[1]-_yCentre,_referencePoint[0]-_xCentre);
  //calculate distance to XYprojected centre of Helix, comparing this with distance to radius around centre gives DistXY
  float DistXY = (_xCentre-xPoint[0])*(_xCentre-xPoint[0]) + (_yCentre-xPoint[1])*(_yCentre-xPoint[1]);
  DistXY = sqrt(DistXY);
  DistXY = fabs(DistXY - _radius);
  
  int nCircles = 0;

  if (fabs(_tanLambda*_radius)>1.0e-20) {
    float xCircles = phi0 - phi -_charge*(xPoint[2]-_referencePoint[2])/(_tanLambda*_radius);
    xCircles = xCircles/_const_2pi;
    int n1,n2;

    if (xCircles >= 0.) {
	n1 = int(xCircles);
	n2 = n1 + 1;
    }
    else {
	n1 = int(xCircles) - 1;
	n2 = n1 + 1;
    }
    
    if (fabs(n1-xCircles) < fabs(n2-xCircles)) {
	nCircles = n1;
    }
    else {
	nCircles = n2;
    }

  }
  
  float DPhi = _const_2pi*((float)nCircles) + phi - phi0;

  zOnHelix = _referencePoint[2] - _charge*_radius*_tanLambda*DPhi;

  float DistZ = fabs(zOnHelix - xPoint[2]);
  float Time;

  if (fabs(_momentum[2]) > 1.0e-20) {
    Time = (zOnHelix - _referencePoint[2])/_momentum[2];
  }
  else {
    Time = _charge*_radius*DPhi/_pxy;
  }

  Distance[0] = DistXY;
  Distance[1] = DistZ;
  Distance[2] = sqrt(DistXY*DistXY+DistZ*DistZ);

  return Time;


}

//When we are not interested in the exact distance, we can check if we are
//already far enough away in XY, before we start calculating in Z as the
//distance will only increase
float HelixClassD::getDistanceToPoint(const std::vector<float>& xPoint, float distCut) {
  //calculate distance to XYprojected centre of Helix, comparing this with distance to radius around centre gives DistXY
  float tempx = xPoint[0]-_xCentre;
  float tempy = xPoint[1]-_yCentre;
  float tempsq = sqrt(tempx*tempx + tempy*tempy);
  float tempdf = tempsq - _radius;
  float DistXY = fabs( tempdf );
  //If this is bigger than distCut, we dont have to know how much bigger this is
  if( DistXY > distCut) {
    return DistXY;
  }

  int nCircles = 0;
  float phi = atan2(tempy,tempx);
  float phi0 = atan2(_referencePoint[1]-_yCentre,_referencePoint[0]-_xCentre);
  float phidiff = phi0-phi;
  float  tempz = xPoint[2] - _referencePoint[2];//Yes referencePoint
  float tanradius = _tanLambda*_radius;
  if (fabs(tanradius)>1.0e-20) {
    float xCircles = (phidiff -_charge*tempz/tanradius)/_const_2pi;
    int n1,n2;
    if (xCircles >= 0.) {
	n1 = int(xCircles);
	n2 = n1 + 1;
    }
    else {
	n1 = int(xCircles) - 1;
	n2 = n1 + 1;
    }
    if (fabs(n1-xCircles) < fabs(n2-xCircles)) {
	nCircles = n1;
    }
    else {
	nCircles = n2;
    }
  }
  float DistZ = - tempz - _charge*tanradius*(_const_2pi*((float)nCircles) - phidiff);
  return sqrt(DistXY*DistXY+DistZ*DistZ);
}//getDistanceToPoint(vector,float)

float HelixClassD::getDistanceToPoint(const float* xPoint, float distCut) {
  std::vector<float> xPosition(xPoint, xPoint + 3 );//We are expecting three coordinates, must be +3, last element is excluded!
  return getDistanceToPoint(xPosition, distCut);
}//getDistanceToPoint(float*,float)



void HelixClassD::setHelixEdges(float * xStart, float * xEnd) {
  for (int i=0; i<3; ++i) {
    _xStart[i] = xStart[i];
    _xEnd[i] = xEnd[i];
  }

}

float HelixClassD::getDistanceToHelix(HelixClassD * helix, float * pos, float * mom) {

  float x01 = getXC();
  float y01 = getYC();
  
  float x02 = helix->getXC();
  float y02 = helix->getYC();
  
  float rad1 = getRadius();
  float rad2 = helix->getRadius();

  float distance = sqrt((x01-x02)*(x01-x02)+(y01-y02)*(y01-y02));

  bool singlePoint = true;

  float phi1 = 0;
  float phi2 = 0;

  if (rad1+rad2<distance) {
    phi1 = atan2(y02-y01,x02-x01);
    phi2 = atan2(y01-y02,x01-x02);
  }
  else if (distance+rad2<rad1) {
    phi1 = atan2(y02-y01,x02-x01);
    phi2 = phi1;
  }
  else if (distance+rad1<rad2) {
    phi1 = atan2(y01-y02,x01-x02);
    phi2 = phi1;
  }
  else {
    singlePoint = false;
    float cosAlpha = 0.5*(distance*distance+rad2*rad2-rad1*rad1)/(distance*rad2);
    float alpha = acos(cosAlpha);
    float phi0 = atan2(y01-y02,x01-x02);
    phi1 = phi0 + alpha;
    phi2 = phi0 - alpha;
  }
  

  float ref1[3];
  float ref2[3];
  for (int i=0;i<3;++i) {
    ref1[i]=_referencePoint[i];
    ref2[i]=helix->getReferencePoint()[i];
  }
  
  float pos1[3];
  float pos2[3];
  float mom1[3];
  float mom2[3];


  if (singlePoint ) {

    float xSect1 = x01 + rad1*cos(phi1);
    float ySect1 = y01 + rad1*sin(phi1);
    float xSect2 = x02 + rad2*cos(phi2);
    float ySect2 = y02 + rad2*sin(phi2);
    float R1 = sqrt(xSect1*xSect1+ySect1*ySect1);
    float R2 = sqrt(xSect2*xSect2+ySect2*ySect2);

    getPointOnCircle(R1,ref1,pos1);
    helix->getPointOnCircle(R2,ref2,pos2);
    
  }
  else {    

    float xSect1 = x02 + rad2*cos(phi1);
    float ySect1 = y02 + rad2*sin(phi1);
    float xSect2 = x02 + rad2*cos(phi2);
    float ySect2 = y02 + rad2*sin(phi2);

//     std::cout << "(xSect1,ySect1)=(" << xSect1 << "," << ySect1 << ")" << std::endl;
//     std::cout << "(xSect2,ySect2)=(" << xSect2 << "," << ySect2 << ")" << std::endl;

    float temp21[3];
    float temp22[3];

    float phiI2  = atan2(ref2[1]-y02,ref2[0]-x02); 
    float phiF21 = atan2(ySect1-y02,xSect1-x02);
    float phiF22 = atan2(ySect2-y02,xSect2-x02);
    float deltaPhi21 = phiF21 - phiI2;
    float deltaPhi22 = phiF22 - phiI2;
    float charge2 = helix->getCharge();
    float pxy2 = helix->getPXY();
    float pz2   = helix->getMomentum()[2];
    if (deltaPhi21 < 0 && charge2 < 0) {
      deltaPhi21 += _const_2pi;
    }
    else if (deltaPhi21 > 0 && charge2 > 0) { 
      deltaPhi21 -= _const_2pi;
    }

    if (deltaPhi22 < 0 && charge2 < 0) {
      deltaPhi22 += _const_2pi;
    }
    else if (deltaPhi22 > 0 && charge2 > 0) { 
      deltaPhi22 -= _const_2pi;
    }

    float time21 = -charge2*deltaPhi21*rad2/pxy2;
    float time22 = -charge2*deltaPhi22*rad2/pxy2;

    float Z21 = ref2[2] + time21*pz2;
    float Z22 = ref2[2] + time22*pz2;

    temp21[0] = xSect1; temp21[1] = ySect1; temp21[2] = Z21;
    temp22[0] = xSect2; temp22[1] = ySect2; temp22[2] = Z22;
    

//     std::cout << "temp21 = " << temp21[0] << " " << temp21[1] << " " << temp21[2] << std::endl;
//     std::cout << "temp22 = " << temp22[0] << " " << temp22[1] << " " << temp22[2] << std::endl;


    float temp11[3];
    float temp12[3];

    float phiI1  = atan2(ref1[1]-y01,ref1[0]-x01); 
    float phiF11 = atan2(ySect1-y01,xSect1-x01);
    float phiF12 = atan2(ySect2-y01,xSect2-x01);
    float deltaPhi11 = phiF11 - phiI1;
    float deltaPhi12 = phiF12 - phiI1;
    float charge1 = _charge;
    float pxy1 = _pxy;
    float pz1   = _momentum[2];
    if (deltaPhi11 < 0 && charge1 < 0) {
      deltaPhi11 += _const_2pi;
    }
    else if (deltaPhi11 > 0 && charge1 > 0) { 
      deltaPhi11 -= _const_2pi;
    }

    if (deltaPhi12 < 0 && charge1 < 0) {
      deltaPhi12 += _const_2pi;
    }
    else if (deltaPhi12 > 0 && charge1 > 0) { 
      deltaPhi12 -= _const_2pi;
    }

    float time11 = -charge1*deltaPhi11*rad1/pxy1;
    float time12 = -charge1*deltaPhi12*rad1/pxy1;

    float Z11 = ref1[2] + time11*pz1;
    float Z12 = ref1[2] + time12*pz1;

    temp11[0] = xSect1; temp11[1] = ySect1; temp11[2] = Z11;
    temp12[0] = xSect2; temp12[1] = ySect2; temp12[2] = Z12;
    
//     std::cout << "temp11 = " << temp11[0] << " " << temp11[1] << " " << temp11[2] << std::endl;
//     std::cout << "temp12 = " << temp12[0] << " " << temp12[1] << " " << temp12[2] << std::endl;

    float Dist1 = 0;
    float Dist2 = 0;

    for (int j=0;j<3;++j) {
      Dist1 += (temp11[j]-temp21[j])*(temp11[j]-temp21[j]);
      Dist2 += (temp12[j]-temp22[j])*(temp12[j]-temp22[j]);
    }

    if (Dist1<Dist2) {
      for (int l=0;l<3;++l) {
	pos1[l] = temp11[l];
	pos2[l] = temp21[l];
      }
    }
    else {
       for (int l=0;l<3;++l) {
	pos1[l] = temp12[l];
	pos2[l] = temp22[l];
      }
    }

  }

  getExtrapolatedMomentum(pos1,mom1);
  helix->getExtrapolatedMomentum(pos2,mom2);

  float helixDistance = 0;

  for (int i=0;i<3;++i) {
    helixDistance += (pos1[i] - pos2[i])*(pos1[i] - pos2[i]);
    pos[i] = 0.5*(pos1[i]+pos2[i]);
    mom[i] = mom1[i]+mom2[i];
  }
  helixDistance = sqrt(helixDistance);

  return helixDistance;

}

void HelixClassD::getExtrapolatedMomentum(float * pos, float * momentum) {

  float phi = atan2(pos[1]-_yCentre,pos[0]-_xCentre);
  if (phi <0.) phi += _const_2pi;
  phi = phi - _phiAtPCA + _phi0;
  momentum[0] = _pxy*cos(phi);
  momentum[1] = _pxy*sin(phi);
  momentum[2] = _momentum[2];


}
