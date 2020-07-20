#ifndef LINECLASS_H
#define LINECLASS_H  
class LineClass {

  public:
  LineClass(float x0,
	    float y0,
	    float z0,
	    float ax,
	    float ay,
	    float az);

  LineClass(float *x0,
	    float *ax);
  
  ~LineClass();

  float * getReferencePoint();
  void setReferencePoint(float *x0);
  float * getDirectionalVector();
  void setDirectionalVector(float *ax);
  float getDistanceToPoint(float * xpoint, float * pos);

 private:

  float _x0[3];
  float _ax[3];


};

#endif
