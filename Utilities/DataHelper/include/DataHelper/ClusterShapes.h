#ifndef ClusterShapes_h
#define ClusterShapes_h


#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <cstdlib>
#include "HelixClass.h"
#include <math.h>


/**
 *    Utility class to derive properties of clusters, such as centre of gravity,
 *    axes of inertia, fits of the cluster shape and so on. All the details are
 *    explained in the documentation of the methods. Several classes of the GSL 
 *    (GNU Scientific Library) are needed in this class.
 *
 *    @authors V. Morgunov (ITEP/DESY), A. Raspereza (DESY), O. Wendt (DESY)
 *    @version $Id: ClusterShapes.h,v 1.14 2007-04-27 13:56:53 owendt Exp $
 *
 */
class ClusterShapes {

public:

  /**
   *    Constructor
   *    @param nhits : number of hits in the cluster
   *    @param a     : amplitudes of elements ('cells') of the cluster. Stored in 
   *                   an array, with one entry for each element ('cell'). Each entry 
   *                   is depending on coordinates x,y,z (Cartesian), which are stored 
   *                   in the arrays x,y,z.
   *    @param x,y,z : array of coordinates corresponding to the array of amplitudes a.
   *
   *
   */
  ClusterShapes(int nhits, float* a, float* x, float* y, float* z);


  /**
   *    Destructor
   */
  ~ClusterShapes();

  /**
   *    Defining errors for Helix Fit
   */
  void setErrors(float *ex, float* ey, float *ez);

  /**
   *   Defining hit types for Helix Fit :
   *   type 1 - cyllindrical detector
   *   type 2 - Z disk detector
   */
  void setHitTypes(int *ityp);


  /**
   * returns the number of elements of the cluster
   */
  int getNumberOfHits();

  /**
   * returns the accumulated amplitude for the whole cluster (just the sum of the 
   * energy in all the entries of the cluster)
   */
  float getTotalAmplitude();

  /**
   * returns an array, which represents a vector from the origin of the
   * coordiante system, i.\ e.\ IP, to the centre of gravity of the cluster. The centre 
   * of gravity is calculated with the energy of the entries of the cluster.
   */
  float* getCentreOfGravity();
  // this is (for now) a pure dummy to allow MarlinPandora development!
  float* getCentreOfGravityErrors();

  /** US spelling of getCentreOfGravity */
  inline float* getCenterOfGravity() { return getCentreOfGravity() ; }
  // this is (for now) a pure dummy to allow MarlinPandora development!
  inline float* getCenterOfGravityErrors() { return getCentreOfGravityErrors() ; }
  
  /**
   * array of the inertias of mass (i.\ e.\ energy) corresponding to the three main axes 
   * of inertia. The array is sorted in ascending order.
   */
  float* getEigenValInertia();
  // this is (for now) a pure dummy to allow MarlinPandora development!
  float* getEigenValInertiaErrors();

  /**
   * array of the three main axes of inertia (9 entries) starting
   * with the axis corresponding to the smallest inertia of mass 
   * (main principal axis). All axes are normalised to a length 
   * of 1.
   */
  float* getEigenVecInertia();
  // this is (for now) a pure dummy to allow MarlinPandora development!
  float* getEigenVecInertiaErrors();

  /**
   * 'mean' width of the cluster perpendicular to the main 
   * principal axis, defined as: 
   * width := sqrt( 1/TotalAmplitude * Sum(a[i]*d[i]*d[i]) ),
   * where d[i] is the distance of the i-th point to the main
   * principal axis.
   */
  float getWidth();

  /**
   * returns the coordinates of the cluster transformed into 
   * the CoG-System.
   * @param xlong  : pointer to an array, where the calculated longitudinal coordiantes
   *                 are stored in.
   * @param xtrans : pointer to an array, where the calculated transversal coordiantes
   *                 are stored in.
   */
  int getEigenSytemCoordinates(float* xlong, float* xtrans);

  /**
   * returns the coordinates and the amplitudes of the cluster
   * transformed into the CoG-System.
   * @param xlong  : pointer to an array, where the calculated longitudinal coordiantes
   *                 are stored in.
   * @param xtrans : pointer to an array, where the calculated transversal coordiantes
   *                 are stored in.
   * @param a      : pointer to an array, where the amplitudes corresponding to the 
   *                 longitudinal and transversal coordiantes are stored in.
   */
  int getEigenSytemCoordinates(float* xlong, float* xtrans, float* a);

  /**
   * performs a least square fit on the shape of an electro-
   * magnetic-shower, which is defined as:
   * A[i] = a * (xl[i]-xl0)^b * exp(-c*(xl[i]-xl0)) * exp(-d*xt[i]),
   * where A[i] is the array of amplitudes, xl[i] is the 
   * coordinate of the actual point along the main principal 
   * axis and xt[i] the coordinate perpendicular to it. The return value 
   * of the method itself is not used at the moment (always returns 0).
   * @param a,b,c,d,xl0  : references to the parameters, which are fitted.
   * @param chi2         : reference to the chi2 of the fit
   * @param xStart       : pointer to the 'initial hit' of the cluster. It is defined 
   *                       as the point with the largest distance to the CoG measured 
   *                       in the direction towards the IP.
   * @param index_xStart : index of the point in the cluster corresponding to xStart
   * @param X0           : radiation length of the detector material. For a composite 
   *                       detector this is meant to be the 'mean' radiation length.
   * @param Rm           : Moliere radius of the the detector material. For a composite 
   *                       detector this is meant to be the 'mean' Moliere radius.
   */
  int fit3DProfile(float& chi2, float& a, float& b, float& c, float& d, float& xl0, 
		   float * xStart, int& index_xStart, float* X0, float* Rm);

  /**
   * returns the chi2 of the fit in the method Fit3DProfile (if simple
   * parametrisation is used)for a given set of parameters a,b,c,d
   * @param a,b,c,d,xl0  : fitted parameters, which have been calculated before
   * @param X0           : radiation length of the detector material. For a composite 
   *                       detector this is meant to be the 'mean' radiation length.
   * @param Rm           : Moliere radius of the the detector material. For a composite 
   *                       detector this is meant to be the 'mean' Moliere radius.
   */
  float getChi2Fit3DProfileSimple(float a, float b, float c, float d, float* X0, 
				  float* Rm);

  /**
   * returns the chi2 of the fit in the method Fit3DProfile (if advanced
   * parametrisation is used) for a given set of parameters E0,a,b,d,t0
   * @param E0,a,b,d,t0 : fitted parameters, which have been calculated before
   * @param X0          : radiation length of the detector material. For a composite 
   *                      detector this is meant to be the 'mean' radiation length.
   * @param Rm          : Moliere radius of the the detector material. For a composite 
   *                      detector this is meant to be the 'mean' Moliere radius.
   */
  float getChi2Fit3DProfileAdvanced(float E0, float a, float b, float d, float t0, 
				    float* X0, float* Rm);

  /**
   * performs a least square fit on a helix path in space, which
   * which is defined as (Cartesian coordiantes):
   *
   * 1. parametrisation:
   * x[i] = x0 + R*cos(b*z[i] + phi0)
   * y[i] = y0 + R*sin(b*z[i] + phi0)
   * z[i] = z[i],
   * where x0,y0,R,b and phi0 are the parameters to be fitted and
   * x[i],y[i],z[i] are the (Cartesian) coordiantes of the space
   * points.
   * 
   * 2. parametrisation:   
   * x[i] = x0 + R*cos(phi)
   * y[i] = y0 + R*sin(phi)
   * z[i] = z0 + b*phi
   * and phi = atan2( y[i]-y0 , x[i]-x0 ),
   * where x0,y0,z0,R and b are the parameters to be fitted and
   * x[i],y[i],z[i] are the (Cartesian) coordiantes of the space
   * points.
   * 
   * The method returns 1 if an error occured and 0 if not.
   *
   * The following output/input parameters are returned/needed:
   *
   * OUTPUTS:
   * @param parameter     : array of parameters to be fitted.
   *                        For parametrisation 1: parameter[5] = {x0,y0,R,b,phi0}
   *                        For parametrisation 2: parameter[5] = {x0,y0,z0,R,b}
   * @param dparameter    : error on the parameters, that means: 
   *                        dparameter[i] = sqrt( CovarMatrix[i][i] )
   * @param chi2          : chi2 of the fit
   * @param distmax       : maximal distance between the points x[i],y[i]
   *                        z[i] and the fitted function
   *
   * INPUTS:
   * @param parametrisation : 1 for first and 2 for second parametrisation (see above)
   * @param max_iter        : maximal number of iterations, before fit cancels
   * @param status_out      : if set to 1, only the initial parameters of
   *                          the fit are calculated and are stored in
   *                          parameter. The entries of dparameter are
   *                          set to 0.0
   */
  int FitHelix(int max_iter, int status_out, int parametrisation,
	       double* parameter, double* dparameter, double& chi2, double& distmax, int direction=1);


  int FitHelix(int max_iter, int status_out, int parametrisation,
	       float* parameter, float* dparameter, float& chi2, float& distmax, int direction=1);

  //here add my functions(variables estimated with detector base)
  //maximum deposit energy of hits
  float getEmax(float* xStart, int& index_xStart, float* X0, float* Rm);

  //shower max of the hits from the shower start hit
  float getsmax(float* xStart, int& index_xStart, float* X0, float* Rm);

  //radius where 90% of the cluster energy exists
  float getxt90(float* xStart, int& index_xStart, float* X0, float* Rm);

  //length where less than 20% of the cluster energy exists
  float getxl20(float* xStart, int& index_xStart, float* X0, float* Rm);

  //for cluster study
  void gethits(float* xStart, int& index_xStart, float* X0, float* Rm, float *okxl, float *okxt, float *oke);
  /**
   * distance to the centre of gravity measured from IP
   * (absolut value of the vector to the centre of gravity)
   */
  inline float radius() { return _radius; }

  /**
   * largest spatial axis length of the ellipsoid derived
   * by the inertia tensor (by their eigenvalues and eigen-
   * vectors)
   */
  inline float getElipsoid_r1() { return _r1; }

  /**
   * medium spatial axis length of the ellipsoid derived
   * by the inertia tensor (by their eigenvalues and eigen-
   * vectors)
   */
  inline float getElipsoid_r2() { return _r2; }

  /**
   * smallest spatial axis length of the ellipsoid derived
   * by the inertia tensor (by their eigenvalues and eigen-   
   * vectors)
   */
  inline float getElipsoid_r3() { return _r3; }

  /**
   * volume of the ellipsoid
   */
  inline float getElipsoid_vol() { return _vol; }

  /**
   * average radius of the ellipsoid (qubic root of volume)
   */
  inline float getElipsoid_r_ave() { return _r_ave; }

  /**
   * density of the ellipsoid defined by: totAmpl/vol
   */
  inline float getElipsoid_density() { return _density; }

  /**
   * eccentricity of the ellipsoid defined by: 
   * Width/r1
   */
  inline float getElipsoid_eccentricity() { return _eccentricity; }

  /**
   * distance from centre of gravity to the point most far 
   * away from IP projected on the main principal axis
   */
  inline float getElipsoid_r_forw() { return _r1_forw; }

  /**
   * distance from centre of gravity to the point nearest 
   * to IP projected on the main principal axis    
   */
  inline float getElipsoid_r_back() { return _r1_back; }

  //Mean of the radius of the hits
  float getRhitMean(float* xStart, int& index_xStart, float* X0, float* Rm);

  //RMS of the radius of the hits
  float getRhitRMS(float* xStart, int& index_xStart, float* X0, float* Rm);



private:

  int _nHits;

  std::vector<float> _aHit;
  std::vector<float> _xHit;
  std::vector<float> _yHit;
  std::vector<float> _zHit;
  std::vector<float> _exHit;
  std::vector<float> _eyHit;
  std::vector<float> _ezHit;
  std::vector<float> _xl;
  std::vector<float> _xt;
  std::vector<float> _t;
  std::vector<float> _s;
  std::vector<int>   _types;

  int   _ifNotGravity=1;
  float _totAmpl=0.0;
  float _radius=0.0;
  float _xgr=0.0;
  float _ygr=0.0;
  float _zgr=0.0;
  float _analogGravity[3]={0.0, 0.0,0.0};

  int   _ifNotWidth=1;
  float _analogWidth=0.0;

  int   _ifNotInertia=1;
  float _ValAnalogInertia[3];
  float _VecAnalogInertia[9];

  int _ifNotEigensystem=1;

  //int   _ifNotElipsoid=1;
  float _r1           =0.0;  // Cluster spatial axis length -- the largest
  float _r2           =0.0;  // Cluster spatial axis length -- less
  float _r3           =0.0;  // Cluster spatial axis length -- less
  float _vol          =0.0;  // Cluster ellipsoid volume
  float _r_ave        =0.0;  // Cluster average radius  (qubic root)
  float _density      =0.0;  // Cluster density
  float _eccentricity =0.0;  // Cluster Eccentricity
  float _r1_forw      =0.0;
  float _r1_back      =0.0;

  void  findElipsoid();
  void  findGravity();
  void  findInertia();
  void  findWidth();
  float findDistance(int i);
  float vecProduct(float * x1, float * x2);
  float vecProject(float * x, float * axis);
  double DistanceHelix(double x, double y, double z, double X0, double Y0, double R0, double bz,
		       double phi0, double * distRPhiZ);
  int transformToEigensystem(float* xStart, int& index_xStart, float* X0, float* Xm);
  float calculateChi2Fit3DProfileSimple(float a, float b, float c, float d);
  float calculateChi2Fit3DProfileAdvanced(float E0, float a, float b, float d, float t0);
  int fit3DProfileSimple(float& chi2, float& a, float& b, float& c, float& d);
  int fit3DProfileAdvanced(float& chi2, double* par_init, double* par, int npar,
			   float* t, float* s, float* E, float E0);

  // private methods for non-linear, multidim. fitting (helix)
  // static int functParametrisation1(const gsl_vector* par, void* data, gsl_vector* f);
  // static int dfunctParametrisation1(const gsl_vector* par, void* d, gsl_matrix* J);
  // static int fdfParametrisation1(const gsl_vector* par, void* d, gsl_vector* f, gsl_matrix* J);


};


#endif
