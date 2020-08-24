#include "DataHelper/LineClass.h"
#include <math.h> 

/*
 * Constructor
 */

LineClass::LineClass(float x0,
		     float y0,
		     float z0,
		     float ax,
		     float ay,
		     float az) {
  _x0[0] = x0;
  _x0[1] = y0;
  _x0[2] = z0;
  _ax[0] = ax;
  _ax[1] = ay;
  _ax[2] = az;
}

LineClass::LineClass(float *x0,
		     float *ax) {

  for (int i=0; i<3; ++i) {
    _x0[i] = x0[i];
    _ax[i] = ax[i];
  }
}

float * LineClass::getReferencePoint() {
  return _x0;
}

float * LineClass::getDirectionalVector() {
  return _ax;
}

void LineClass::setReferencePoint(float *x0) {
  for (int i=0; i<3; ++i)
    _x0[i] = x0[i];
}
 
void LineClass::setDirectionalVector(float *ax) {
  for (int i=0; i<3; ++i)
    _ax[i] = ax[i];
 
}

float LineClass::getDistanceToPoint(float * xpoint, float * pos) {
  
  float dif[3];
  float prod = 0;
  float den = 0;
  for (int i=0; i<3; ++i) {
    dif[i] = xpoint[i] - _x0[i];
    prod += _ax[i]*dif[i];
    den += _ax[i]*_ax[i];
  }
  float time = prod/fmax(1e-10,den);

  float dist = 0.0;
  for (int i=0; i<3; ++i) {
    pos[i] = _x0[i] + _ax[i]*time;
    dist += (xpoint[i]-pos[i])*(xpoint[i]-pos[i]);
  }
  dist = sqrt(dist);

  return dist;



}
