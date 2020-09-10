

/* -*- Mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */

#include "DataHelper/ClusterShapes.h"

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
//#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_integration.h>
//#include <gsl/gsl_rng.h>
//#include <gsl/gsl_sf_pow_int.h>


// #################################################
// #####                                       #####
// #####  Additional Structures and Functions  #####
// #####                                       #####
// #################################################
//=============================================================================

struct data {
  int n;
  float* x;
  float* y;
  float* z;
  float* ex;
  float* ey;
  float* ez;
};

//=============================================================================

// Gammafunction
double G(double x) {

  return gsl_sf_gamma(x);
  
}

//=============================================================================

// inverse Gammafunction
double invG(double x) {
  
  return gsl_sf_gammainv(x);
    
}

//=============================================================================

// Integral needed for deriving the derivative of the Gammafunction
double Integral_G(double x, void* params) {
  double a = *(double*)params;
  double f = exp(-x) * pow(x,a-1) * log(x);
  return f;
}
    
//=============================================================================

double DinvG(double x) {

  int workspace_size = 1000;
  double abs_error = 0;
  double rel_error = 1e-6;
  double result = 0.0;
  double error = 0.0;
  
  gsl_integration_workspace* w  = gsl_integration_workspace_alloc(workspace_size);
  gsl_function F;
  F.function = &Integral_G;
  F.params = &x;
    
  /*int status=*/gsl_integration_qagiu(&F,0,abs_error,rel_error,workspace_size,w,
				 &result,&error); 

  // debug  
  /*  
      printf ("Numeric Integration : \n");
      printf ("parameter of integration = % .18f\n", x);
      printf ("status of integration    = %d \n"   , status);
      printf ("result                   = % .18f\n", result);
      printf ("estimated error          = % .18f\n", error);
      printf ("intervals                =  %d\n\n", w->size);
  */
  
  double G2 = pow(gsl_sf_gamma(x),2); 
  double DG = result;
  
  gsl_integration_workspace_free(w);
  
  return -DG/G2;

}

//=============================================================================

int ShapeFitFunct(const gsl_vector* par, void* d, gsl_vector* f) {

  // Used for shape fitting. Function to fit: 
  //
  // a[i](t[i],s[i]) =
  // 
  // E0 * b * 1/Gamma(a) * ( b * (t[i] - t0) )^(a-1) * exp(-b*(t[i] - t0)) * exp(-d*s[i])
  //
  // Function to minimise:
  //
  // f0[i] =  E0 * b * 1/Gamma(a) * 
  //        ( b * (t[i] - t0) )^(a-1) * exp(-b*(t[i] - t0)) * exp(-d*s[i]) - a[i]
  //


  //  float E0   = gsl_vector_get(par,0);
  float A    = gsl_vector_get(par,0);
  float B    = gsl_vector_get(par,1);
  float D    = gsl_vector_get(par,2);
  float t0   = gsl_vector_get(par,3);
  int n      = ((struct data*)d)->n;
  float* t   = ((struct data*)d)->x;
  float* s   = ((struct data*)d)->y;
  float* a   = ((struct data*)d)->z; // amplitude stored in z[i]
  float fi   = 0.0;


  for (int i(0); i < n; i++) {
    fi = /*E0 * */ B * invG(A) * pow(B*(t[i]-t0),A-1) * exp(-B*(t[i]-t0)) * exp(-D*s[i]) 
      - a[i];
    gsl_vector_set(f,i,fi);
  }

  return GSL_SUCCESS;
}

//=============================================================================

int dShapeFitFunct(const gsl_vector* par, void* d, gsl_matrix* J) {

  // Used for shape fitting
  //float E0  = gsl_vector_get(par,0);
  float A   = gsl_vector_get(par,0);
  float B   = gsl_vector_get(par,1);
  float D   = gsl_vector_get(par,2);
  float t0  = gsl_vector_get(par,3);
  int n     = ((struct data*)d)->n;
  float* t  = ((struct data*)d)->x;
  float* s  = ((struct data*)d)->y;


  // calculate Jacobi's matrix J[i][j] = dfi/dparj, but here only one dimension

  for (int i(0); i < n; i++) {

    /*    
	  gsl_matrix_set(J,i,0,B * invG(A) * pow(B*(t[i]-t0),A-1) * exp(-B*(t[i]-t0)) 
	  * exp(-D*s[i]) );
	  */

    gsl_matrix_set(J,i,0,( /* E0 * */ B * invG(A) * log(B*(t[i]-t0))*pow(B*(t[i]-t0),A-1) * 
			   exp(-B*(t[i]-t0)) + DinvG(A) * /* E0 * */ B * pow(B*(t[i]-t0),A-1) *
			   exp(-B*(t[i]-t0))
			   ) * exp(-D*s[i]));
		         
    gsl_matrix_set(J,i,1,( /* E0 * */ invG(A) * pow(B*(t[i]-t0),A-1) * exp(-B*(t[i]-t0)) +
			   /* E0 * */ invG(A) * (A-1) * B * (t[i]-t0) * pow(B*(t[i]-t0),A-2) * 
			   exp(-B*(t[i]-t0)) -
			   /* E0 * */ B * invG(A) * (t[i]-t0) * pow(B*(t[i]-t0),A-1) * 
			   exp(-B*(t[i]-t0)) 
			   ) * exp(-D*s[i]));
                         
    gsl_matrix_set(J,i,2,-/* E0 * */ B * invG(A) * s[i] * pow(B*(t[i]-t0),A-1) *
		   exp(-B*(t[i]-t0)) * exp(-D*s[i]));
       		         
    gsl_matrix_set(J,i,3,(-/* E0 * */ pow(B,2) * invG(A) * (A-1) * pow(B*(t[i]-t0),A-2) * 
			  exp(-B*(t[i]-t0)) +
			  /* E0 * */ pow(B,2) * invG(A) * pow(B*(t[i]-t0),A-1) * 
			  exp(-B*(t[i]-t0))
			  ) * exp(-D*s[i]));
 		         

  }
  
  return GSL_SUCCESS;
}

//=============================================================================

int fdfShapeFitFunct(const gsl_vector* par, void* d, gsl_vector* f, gsl_matrix* J) {

  //     For helix fitting
  ShapeFitFunct(par, d, f);
  dShapeFitFunct(par, d, J);

  return GSL_SUCCESS;

}

//=============================================================================

int signum(float x) {

  // computes the signum of x. Needed for the 3. parametrisation

  if ( x >= 0 ) return 1; // x == 0 is taken as positive
  else return -1;

}

//=============================================================================

int functParametrisation1(const gsl_vector* par, void* d, gsl_vector* f) {

  //     For helix fitting
  // calculate fit function f0[i] = 
  // ( (x0 + R*cos(b*z[i] + phi0)) - x[i] ) for i = 0 to n-1
  //                    and f1[i] = 
  // ( (y0 + R*sin(b*z[i] + phi0)) - y[i] ) for i = n to dim*n - 1
  // That means, minimise the two functions f0[i] and f1[i]

  float x0   = gsl_vector_get(par,0);
  float y0   = gsl_vector_get(par,1);
  float R    = gsl_vector_get(par,2);
  float b    = gsl_vector_get(par,3);
  float phi0 = gsl_vector_get(par,4);
  int n    = ((struct data*)d)->n;
  float* x = ((struct data*)d)->x;
  float* y = ((struct data*)d)->y;
  float* z = ((struct data*)d)->z;
  //float* ex = ((struct data*)d)->ex;
  //float* ey = ((struct data*)d)->ey;
  //float* ez = ((struct data*)d)->ez;
  float fi = 0.0;

  // first dimension
  for (int i(0); i < n; i++) {
    fi = (x0 + R*cos(b*z[i] + phi0)) - x[i];
    //    float err = sqrt(ex[i]*ex[i]+R*b*R*b*sin(b*z[i] + phi0)*R*b*R*b*sin(b*z[i] + phi0)*ez[i]*ez[i]);
    //    fi = fi/err;
    gsl_vector_set(f,i,fi);
  }
  // second dimension
  for (int i(0); i < n; i++) {
    fi = (y0 + R*sin(b*z[i] + phi0)) - y[i];
    //    float err = sqrt(ey[i]*ey[i]+R*b*R*b*cos(b*z[i] + phi0)*R*b*R*b*cos(b*z[i] + phi0)*ez[i]*ez[i]);
    //    fi = fi/err;  
    gsl_vector_set(f,i+n,fi);
  }

  return GSL_SUCCESS;
}

//=============================================================================

int dfunctParametrisation1(const gsl_vector* par, void* d, gsl_matrix* J) {

  //     For helix fitting
  float R    = gsl_vector_get(par,2);
  float b    = gsl_vector_get(par,3);
  float phi0 = gsl_vector_get(par,4);

  int n    = ((struct data*)d)->n;
  float* z = ((struct data*)d)->z;
  //float* ex = ((struct data*)d)->ex;
  //float* ey = ((struct data*)d)->ey;
  //float* ez = ((struct data*)d)->ez;


  // calculate Jacobi's matrix J[i][j] = dfi/dparj

  // part of Jacobi's matrix corresponding to first dimension
  for (int i(0); i < n; i++) {
    //    float err = sqrt(ex[i]*ex[i]+R*b*R*b*sin(b*z[i] + phi0)*R*b*R*b*sin(b*z[i] + phi0)*ez[i]*ez[i]);
    gsl_matrix_set(J,i,0,1);
    gsl_matrix_set(J,i,1,0);
    gsl_matrix_set(J,i,2,cos(b*z[i]+phi0));
    gsl_matrix_set(J,i,3,-z[i]*R*sin(b*z[i]+phi0));
    gsl_matrix_set(J,i,4,-R*sin(b*z[i]+phi0));

  }

  // part of Jacobi's matrix corresponding to second dimension
  for (int i(0); i < n; i++) {
    //    float err = sqrt(ey[i]*ey[i]+R*b*R*b*cos(b*z[i] + phi0)*R*b*R*b*cos(b*z[i] + phi0)*ez[i]*ez[i]);
    gsl_matrix_set(J,i+n,0,0);
    gsl_matrix_set(J,i+n,1,1);
    gsl_matrix_set(J,i+n,2,sin(b*z[i]+phi0));
    gsl_matrix_set(J,i+n,3,z[i]*R*cos(b*z[i]+phi0));
    gsl_matrix_set(J,i+n,4,R*cos(b*z[i]+phi0));
    
  }
  
  return GSL_SUCCESS;
}

//=============================================================================

int fdfParametrisation1(const gsl_vector* par, void* d, gsl_vector* f, gsl_matrix* J) {

  //     For helix fitting
  functParametrisation1(par, d, f);
  dfunctParametrisation1(par, d, J);

  return GSL_SUCCESS;

}

//=============================================================================

int functParametrisation2(const gsl_vector* par, void* d, gsl_vector* f) {

  //     For helix fitting
  // calculate fit function f0[i] = 
  // ( (x0 + R*cos(phi)) - x[i] ) for i = 0 to n-1
  //                        f1[i] = 
  // ( (y0 + R*sin(phi)) - y[i] ) for i = n to dim*n - 1
  //                    and f2[i] =
  // ( (z0 + b*phi     ) - z[i] )
  // That means, minimise the three functions f0[i], f1[i] and f2[i]

  float x0   = gsl_vector_get(par,0);
  float y0   = gsl_vector_get(par,1);
  float z0   = gsl_vector_get(par,2);
  float R    = gsl_vector_get(par,3);
  float b    = gsl_vector_get(par,4);
  int n    = ((struct data*)d)->n;
  float* x = ((struct data*)d)->x;
  float* y = ((struct data*)d)->y;
  float* z = ((struct data*)d)->z;
  float fi = 0.0;
  float phii = 0.0;

  // first dimension
  for (int i(0); i < n; i++) {
    phii = atan2( y[i]-y0, x[i]-x0 );
    fi = (x0 + R*cos(phii)) - x[i];
    gsl_vector_set(f,i,fi);
  }
  // second dimension
  for (int i(0); i < n; i++) {
    phii = atan2( y[i]-y0, x[i]-x0 );
    fi = (y0 + R*sin(phii)) - y[i];  
    gsl_vector_set(f,i+n,fi);
  }
  // third dimension
  for (int i(0); i < n; i++) {
    phii = atan2( y[i]-y0, x[i]-x0 );
    fi = (z0 + b*phii) - z[i];  
    gsl_vector_set(f,i+2*n,fi);
  }

  return GSL_SUCCESS;
}

//=============================================================================

int dfunctParametrisation2(const gsl_vector* par, void* d, gsl_matrix* J) {

  //     For helix fitting
  float x0   = gsl_vector_get(par,0);
  float y0   = gsl_vector_get(par,1);
  float R    = gsl_vector_get(par,3);
  float b    = gsl_vector_get(par,4);
  int n    = ((struct data*)d)->n;
  float* x = ((struct data*)d)->x;
  float* y = ((struct data*)d)->y;
  float phii = 0.0;

  // calculate Jacobi's matrix J[i][j] = dfi/dparj

  // part of Jacobi's matrix corresponding to first dimension
  for (int i(0); i < n; i++) {
    
    phii = atan2( y[i]-y0, x[i]-x0 );

    gsl_matrix_set(J,i,0,1 - R*sin(phii)*
		   ( (y[i]-y0)/( (x[i]-x0)*(x[i]-x0) + (y[i]-y0)*(y[i]-y0) ) ) );
    gsl_matrix_set(J,i,1,R*sin(phii)*
		   ( (x[i]-x0)/( (x[i]-x0)*(x[i]-x0) + (y[i]-y0)*(y[i]-y0) ) ) );
    gsl_matrix_set(J,i,2,0);
    gsl_matrix_set(J,i,3,cos(phii));
    gsl_matrix_set(J,i,4,0);

  }
  
  // part of Jacobi's matrix corresponding to second dimension
  for (int i(0); i < n; i++) {

    phii = atan2( y[i]-y0, x[i]-x0 );
    
    gsl_matrix_set(J,i+n,0,R*cos(phii)*
		   ( (y[i]-y0)/( (x[i]-x0)*(x[i]-x0) + (y[i]-y0)*(y[i]-y0) ) ) );
    gsl_matrix_set(J,i+n,1,1 + R*cos(phii)*
		   ( (x[i]-x0)/( (x[i]-x0)*(x[i]-x0) + (y[i]-y0)*(y[i]-y0) ) ) );
    gsl_matrix_set(J,i+n,2,0);
    gsl_matrix_set(J,i+n,3,sin(phii));
    gsl_matrix_set(J,i+n,4,0);
    
  }

  // part of Jacobi's matrix corresponding to third dimension
  for (int i(0); i < n; i++) {

    phii = atan2( y[i]-y0, x[i]-x0 );
    
    gsl_matrix_set(J,i+2*n,0,b*
		   ( (y[i]-y0)/( (x[i]-x0)*(x[i]-x0) + (y[i]-y0)*(y[i]-y0) ) ) );
    gsl_matrix_set(J,i+2*n,1,b*
		   ( (x[i]-x0)/( (x[i]-x0)*(x[i]-x0) + (y[i]-y0)*(y[i]-y0) ) ) );
    gsl_matrix_set(J,i+2*n,2,1);
    gsl_matrix_set(J,i+2*n,3,0);
    gsl_matrix_set(J,i+2*n,4,phii);
    
  }
  
  return GSL_SUCCESS;
}

//=============================================================================

int fdfParametrisation2(const gsl_vector* par, void* d, gsl_vector* f, gsl_matrix* J) {

  //     For helix fitting
  functParametrisation2(par, d, f);
  dfunctParametrisation2(par, d, J);

  return GSL_SUCCESS;

}

//=============================================================================

int functParametrisation3(const gsl_vector* par, void* d, gsl_vector* f) {

  //     For helix fitting
  // calculate fit function f0[i] = 
  // ( ( ( (1/omega) - d0 )*sin(Phi0) + ( 1/fabs(omega) )*cos( ( -omega/sqrt(1+tanL^2) )*s + Phi0 +( (omega*pi)/(2*fabs(omega)) ) ) ) - x[i] ) for i = 0 to n-1
  //                        f1[i] = 
  // ( ( (-1.0)*( (1/omega) - d0 )*cos(Phi0) + ( 1/fabs(omega) )*sin( ( -omega/sqrt(1+tanL^2) )*s + Phi0 +( (omega*pi)/(2*fabs(omega)) ) ) ) - y[i] ) for i = n to dim*n - 1
  //                    and f2[i] =
  // ( ( z0 + (tanL/sqrt(1+tanL^2))*s ) - z[i] )
  // That means, minimise the three functions f0[i], f1[i] and f2[i]

  double z0    = gsl_vector_get(par,0);
  double Phi0  = gsl_vector_get(par,1);
  double omega = gsl_vector_get(par,2);
  double d0    = gsl_vector_get(par,3);
  double tanL  = gsl_vector_get(par,4);
  int n    = ((struct data*)d)->n;
  float* x = ((struct data*)d)->x;
  float* y = ((struct data*)d)->y;
  float* z = ((struct data*)d)->z;
  double phii = 0.0;
  double fi = 0.0;
  double si = 0.0;

  // first dimension
  for (int i(0); i < n; i++) {
    phii = atan2( ( ((double)y[i]) + ((1/omega) - d0 )*cos(Phi0) ), ( ((double)x[i]) - ((1/omega) - d0 )*sin(Phi0) ) );
    fi = ( ( (1/omega) - d0 )*sin(Phi0) + ( 1/fabs(omega) )*cos(phii) ) - ((double)x[i]);
    gsl_vector_set(f,i,fi);
  }
  // second dimension
  for (int i(0); i < n; i++) {
    phii = atan2( ( ((double)y[i]) + ((1/omega) - d0 )*cos(Phi0) ), ( ((double)x[i]) - ((1/omega) - d0 )*sin(Phi0) ) );
    fi = ( (-1.0)*( (1/omega) - d0 )*cos(Phi0) + ( 1/fabs(omega) )*sin(phii) ) - ((double)y[i]);
    gsl_vector_set(f,i+n,fi);
  }
  // third dimension
  for (int i(0); i < n; i++) {
    phii = atan2( ( ((double)y[i]) + ((1/omega) - d0 )*cos(Phi0) ), ( ((double)x[i]) - ((1/omega) - d0 )*sin(Phi0) ) );
    si = (-1.0)*( (sqrt(1+pow(tanL,2)))/omega )*(phii - Phi0 - (omega*M_PI)/(2*fabs(omega)));
    fi = ( z0 + (tanL/sqrt(1+pow(tanL,2)))*si ) - ((double)z[i]);
    gsl_vector_set(f,i+2*n,fi);
  }
  
  return GSL_SUCCESS;

}

//=============================================================================

int dfunctParametrisation3(const gsl_vector* par, void* d, gsl_matrix* J) {

 
  //     For helix fitting
  // double z0    = gsl_vector_get(par,0); // not needed
  double Phi0  = gsl_vector_get(par,1);
  double omega = gsl_vector_get(par,2);
  double d0    = gsl_vector_get(par,3);
  double tanL  = gsl_vector_get(par,4);
  int n    = ((struct data*)d)->n;
  float* x = ((struct data*)d)->x;
  float* y = ((struct data*)d)->y;
  // float* z = ((struct data*)d)->z; // not needed
  double phii = 0.0;
  double si = 0.0;

  // calculate Jacobi's matrix J[i][j] = dfi/dparj

  // part of Jacobi's matrix corresponding to first dimension
  for (int i(0); i < n; i++) {
    phii = atan2( ( ((double)y[i]) + ((1/omega) - d0 ) * cos(Phi0) ), ( ((double)x[i]) - ((1/omega) - d0 )*sin(Phi0) ) );
    si = (-1.0)*( (sqrt(1+pow(tanL,2)))/omega )*(phii - Phi0 - (omega*M_PI)/(2*fabs(omega)));

    gsl_matrix_set(J,i,0,0);
    gsl_matrix_set(J,i,1,((1/omega) - d0) * cos(Phi0) - (1/fabs(omega)) * sin(phii) );
    gsl_matrix_set(J,i,2,((-1.0)*sin(Phi0))/pow(omega,2) - ( (signum(omega))/(pow(fabs(omega),2)) ) * cos(phii) - 
		   (1/fabs(omega)) * sin(phii) * ( ( ((-1.0)/sqrt(1+pow(tanL,2)))*si) + (M_PI)/(2*fabs(omega)) - (signum(omega)*omega*M_PI)/(2*pow(fabs(omega),2)) ) );
    gsl_matrix_set(J,i,3,(-1.0)*sin(Phi0));
    gsl_matrix_set(J,i,4,((-1.0)/fabs(omega))*sin(phii) * ( (omega*tanL*si)/sqrt(pow(1+pow(tanL,2),3)) ) );

  }
  
  // part of Jacobi's matrix corresponding to second dimension
  for (int i(0); i < n; i++) {
    phii = atan2( ( ((double)y[i]) + ((1/omega) - d0 )*cos(Phi0) ), ( ((double)x[i]) - ((1/omega) - d0 )*sin(Phi0) ) );
    si = (-1.0)*( (sqrt(1+pow(tanL,2)))/omega )*(phii - Phi0 - (omega*M_PI)/(2*fabs(omega)));

    gsl_matrix_set(J,i+n,0,0);
    gsl_matrix_set(J,i+n,1,((1/omega) - d0)*sin(Phi0) + (1/fabs(omega)) * cos(phii) );
    gsl_matrix_set(J,i+n,2,cos(Phi0)/pow(omega,2) + ( (signum(omega))/(pow(fabs(omega),2)) ) * sin(phii) + 
		   (1/fabs(omega))*cos(phii) * ( ( ((-1.0)/sqrt(1+pow(tanL,2)))*si) + (M_PI)/(2*fabs(omega)) - (signum(omega)*omega*M_PI)/(2*pow(fabs(omega),2)) ) );
    gsl_matrix_set(J,i+n,3,cos(Phi0));
    gsl_matrix_set(J,i+n,4,(1/fabs(omega))*cos(phii) * ( (omega*tanL*si)/sqrt(pow(1+pow(tanL,2),3)) ) );
    
  }

  // part of Jacobi's matrix corresponding to third dimension
  for (int i(0); i < n; i++) {
    phii = atan2( ( ((double)y[i]) + ((1/omega) - d0 )*cos(Phi0) ), ( ((double)x[i]) - ((1/omega) - d0 )*sin(Phi0) ) );
    si = (-1.0)*( (sqrt(1+pow(tanL,2)))/omega )*(phii - Phi0 - (omega*M_PI)/(2*fabs(omega)));

    gsl_matrix_set(J,i+2*n,0,1.0);
    gsl_matrix_set(J,i+2*n,1,0);
    gsl_matrix_set(J,i+2*n,2,0);
    gsl_matrix_set(J,i+2*n,3,0);
    gsl_matrix_set(J,i+2*n,4,si/sqrt(1+pow(tanL,2)) - (pow(tanL,2)*si)/sqrt(pow(1+pow(tanL,2),3)) );
    
  }
  
  return GSL_SUCCESS;

}

//=============================================================================

int fdfParametrisation3(const gsl_vector* par, void* d, gsl_vector* f, gsl_matrix* J) {
  
  //     For helix fitting
  functParametrisation3(par, d, f);
  dfunctParametrisation3(par, d, J);
  
  return GSL_SUCCESS;

}

//=============================================================================




// ##########################################
// #####                                #####
// #####   Constructor and Destructor   #####
// #####                                #####
// ##########################################

//=============================================================================

ClusterShapes::ClusterShapes(int nhits, float* a, float* x, float* y, float* z):

  _nHits(nhits),
  _aHit  (nhits, 0.0),
  _xHit  (nhits, 0.0),
  _yHit  (nhits, 0.0),
  _zHit  (nhits, 0.0),
  _exHit (nhits, 1.0),
  _eyHit (nhits, 1.0),
  _ezHit (nhits, 1.0),
  _xl    (nhits, 0.0),
  _xt    (nhits, 0.0),
  _t     (nhits, 0.0),
  _s     (nhits, 0.0),
  _types(nhits, 1), // all hits are assumed to be "cylindrical"

  _ifNotGravity    (1),
  _ifNotWidth      (1),
  _ifNotInertia    (1),
  _ifNotEigensystem(1)
{

  for (int i(0); i < nhits; ++i) {
    _aHit[i] = a[i];
    _xHit[i] = x[i];
    _yHit[i] = y[i];
    _zHit[i] = z[i];
  }


}


//=============================================================================

ClusterShapes::~ClusterShapes() {

}

//=============================================================================

void ClusterShapes::setErrors(float *ex, float *ey, float *ez) {

  for (int i=0; i<_nHits; ++i) {
    _exHit[i] = ex[i];
    _eyHit[i] = ey[i];
    _ezHit[i] = ez[i];
  }	

}

void ClusterShapes::setHitTypes(int* typ) {
  for (int i=0; i<_nHits; ++i) {
    _types[i] = typ[i];

  }
  
}



// ##########################################
// #####                                #####
// #####        public methods          #####
// #####                                #####
// ##########################################

//=============================================================================

int ClusterShapes::getNumberOfHits() {
  return _nHits;
}

//=============================================================================

float ClusterShapes::getTotalAmplitude() {
  if (_ifNotGravity == 1) findGravity();
  return _totAmpl;
}

//=============================================================================

float* ClusterShapes::getCentreOfGravity() {
  if (_ifNotGravity == 1) findGravity() ;
  return &_analogGravity[0] ;
}
float* ClusterShapes::getCentreOfGravityErrors() {
  // this is a pure dummy to allow MarlinPandora development!
  if (_ifNotGravity == 1) findGravity() ;
  return &_analogGravity[0] ;
}

//=============================================================================

float* ClusterShapes::getEigenValInertia() {
  if (_ifNotInertia == 1) findInertia();
  return &_ValAnalogInertia[0] ;
}
float* ClusterShapes::getEigenValInertiaErrors() {
  // this is a pure dummy to allow MarlinPandora development!
  if (_ifNotInertia == 1) findInertia();
  return &_ValAnalogInertia[0] ;
}

//=============================================================================

float* ClusterShapes::getEigenVecInertia() {
  if (_ifNotInertia == 1) findInertia();
  return &_VecAnalogInertia[0] ;
}
float* ClusterShapes::getEigenVecInertiaErrors() {
  // this is a pure dummy to allow MarlinPandora development!
  if (_ifNotInertia == 1) findInertia();
  return &_VecAnalogInertia[0] ;
}

//=============================================================================

float ClusterShapes::getWidth() {
  if (_ifNotWidth == 1) findWidth();
  return _analogWidth;
}

//=============================================================================

int ClusterShapes::getEigenSytemCoordinates(float* xlong, float* xtrans) {

  float xStart[3];
  int index_xStart;


  // NOT SAVE, change to class variables !!!!!
  float X0[2]={3.50,17.57};   //in mm. //this is the exact value of tungsten and iron
  float Rm[2]={9.00,17.19};   //in mm. need to change to estimate correctly 

  if (_ifNotEigensystem == 1) transformToEigensystem(xStart,index_xStart,X0,Rm);

  for (int i = 0; i < _nHits; ++i) {
    xlong[i]  = _xl[i];
    xtrans[i] = _xt[i];
  }

  return 0; // no error messages at the moment

}

//=============================================================================

int ClusterShapes::getEigenSytemCoordinates(float* xlong, float* xtrans, float* a) {

  float xStart[3];
  int index_xStart;

  // NOT SAVE, change to class variables !!!!!
  float X0[2]={3.50,17.57};   //in mm. //this is the exact value of tungsten and iron
  float Rm[2]={9.00,17.19};   //in mm. need to change to estimate correctly 

  if (_ifNotEigensystem == 1) transformToEigensystem(xStart,index_xStart,X0,Rm);

  for (int i = 0; i < _nHits; ++i) {
    xlong[i]  = _xl[i];
    xtrans[i] = _xt[i];
    a[i] = _aHit[i];
  }

  return 0; // no error messages at the moment

}

//=============================================================================

int ClusterShapes::fit3DProfile(float& chi2, float& E0, float& A, float& B, float& D,
				float& xl0, float* xStart, int& index_xStart,
				float* X0, float* Rm) {

  const int npar = 4;

  if (_ifNotEigensystem == 1){
    transformToEigensystem(xStart,index_xStart,X0,Rm);
  }
  
  float* E = new float[_nHits];

  //doesn't fit when _nHits==0
  //std::cout << "_nhits " << _nHits << std::endl;
  if(_nHits-1 < npar){  //can't fit because number of degrees of freedom is small
    E0=0.0;

    delete[] E;
    int result = 0;  // no error handling at the moment
    return result;
  }

  double par_init[npar];
  for (int i = 0; i < npar; ++i) par_init[i] = 0.0; // initialise

  float E0_init = 0.0;
  float A_init  = 0.0;
  float B_init  = 0.5; // empirically
  float D_init  = 0.06;    //0.5; //if Rm includes 90% of shower energy, absorption length has 63% of shower energy
  float t0_init = -10.0/X0[0]; // shift xl0_ini is assumed to be -10.0 without a reason ????
  
  float E0_tmp = 0.0;
  //  int i_max = 0;
  //float t_max = 0.0;
  for (int i = 0; i < _nHits; ++i) {
    if (E0_tmp < _aHit[i]) {
      E0_tmp = _aHit[i]; 
      //i_max = i;
    }
    E0_init += _aHit[i];
  }

  // first definition
  //t_max = _xl[i_max]/X0;
  //A_init = t_max*B_init + 1.0;

  // second definition
  //t_max = (1.0/3.0) * (t[_nHits-1] - t[0]);
  //A_init = t_max*B_init + 1.0;

  // third definition
  float Ec = X0[0] * 0.021/Rm[0];
  A_init =  B_init * log(E0_init/Ec) + 0.5 * B_init + 1.0; // (+0.5 for Photons initiated showers)

  // par_init[0] = E0_init;
  E0 = E0_init;
  for (int i = 0; i < _nHits; ++i) E[i] = _aHit[i]/E0_init;
  par_init[0] = A_init; // 2.0
  par_init[1] = B_init;
  par_init[2] = D_init;
  par_init[3] = t0_init;
 
  // debug

  // std::cout << "E0_init : " <<  E0_init << "\t" << "A_init : " << A_init << "\t" 
  //     << "B_init : " <<  B_init << "\t" << "D_init : " << D_init << "\t" 
  //     << "xl0_init : "  << t0_init*X0 << "\t" << "X0 : " << X0 
  //     << "\t" << "t_max : " << t_max << std::endl << std::endl;


  double t0 = xl0/X0[0];   //probably t0 is in ecal
  double par[npar];
  par[0] = A;
  par[1] = B;
  par[2] = D;
  par[3] = t0;

  fit3DProfileAdvanced(chi2,par_init,par,npar,&_t[0],&_s[0],E,E0);

  A   = par[0];
  B   = par[1];
  D   = par[2];
  xl0 = par[3] * X0[0];   //probably shower start is in ecal

  delete[] E;

  int result = 0;  // no error handling at the moment
  return result;
}

//=============================================================================

float ClusterShapes::getChi2Fit3DProfileSimple(float a, float b, float c, float d,
					       float* X0, float* Rm) {

  float chi2 = 0.0;

  float xStart[3];
  int index_xStart;

  if (_ifNotEigensystem == 1) transformToEigensystem(xStart,index_xStart,X0,Rm);

  chi2 = calculateChi2Fit3DProfileSimple(a,b,c,d);
  
  return chi2;
  
}

//=============================================================================

float ClusterShapes::getChi2Fit3DProfileAdvanced(float E0, float a, float b, float d,
						 float t0, float* X0, float* Rm) {

  float chi2 = 0.0;

  float xStart[3];
  int index_xStart;

  if (_ifNotEigensystem == 1) transformToEigensystem(xStart,index_xStart,X0,Rm);

  chi2 = calculateChi2Fit3DProfileAdvanced(E0,a,b,d,t0);

  return chi2;

}

//=============================================================================

int ClusterShapes::FitHelix(int max_iter, int status_out, int parametrisation,
			    float* parameter, float* dparameter, float& chi2, 
			    float& distmax, int direction) {


  // Modified by Hengne Li @ LAL
  
  double parameterdb[5];
  double dparameterdb[5];
  double chi2db;
  double distmaxdb;
  for ( int i=0; i<5; i++ ){
    parameterdb[i] = double(parameter[i]);
    dparameterdb[i] = double(dparameter[i]);
  }
  chi2db = double(chi2);
  distmaxdb = double(distmax);
  
  int returnvalue = FitHelix(max_iter,status_out,parametrisation,parameterdb,dparameterdb,chi2db,distmaxdb,direction);
  
  for ( int i=0; i<5; i++ ){
    parameter[i] = float(parameterdb[i]);
    dparameter[i] = float(dparameterdb[i]);
  }
  chi2 = float(chi2db);
  distmax = float(distmaxdb);
  
  return returnvalue ;

}

//=============================================================================

int ClusterShapes::FitHelix(int max_iter, int status_out, int parametrisation,
			    double* parameter, double* dparameter, double& chi2, 
			    double& distmax, int direction) {

  // FIXME: version with double typed parameters needed 2006/06/10 OW
  
  if (_nHits < 3) {
    std::cout << "ClusterShapes : helix fit impossible, two few points" ;
    std::cout << std::endl;
    for (int i = 0; i < 5; ++i) {
      parameter[i] = 0.0;
      dparameter[i] = 0.0;
    }
    return 1;
  }

  // find initial parameters

  double Rmin = 1.0e+10;
  double Rmax = -1.0;
  int i1 = -1;

  // 1st loop  
  for (int i(0); i < _nHits; ++i) {
    double Rz = sqrt(_xHit[i]*_xHit[i] + _yHit[i]*_yHit[i]);
    if (Rz < Rmin) {
      Rmin = Rz;
      i1 = i;
    }
    if (Rz > Rmax) {
      Rmax = Rz;
    }

  }



  // debug
  /*
    for (int i(0); i < _nHits; ++i) std::cout << i << "  " << _xHit[i] << "  " << _yHit[i] << "  " << _zHit[i] << std::endl;
    std::cout << std::endl << Rmin << "  " <<  Rmax << "  " << i1 << std::endl;
  */





  // 2nd loop
  double Upper = Rmin + 1.1*(Rmax-Rmin);
  double Lower = Rmin + 0.9*(Rmax-Rmin);
  double dZmin  = 1.0e+20;

  int i3 = -1 ;

  for (int i(0); i < _nHits; ++i) {
    double Rz = sqrt(_xHit[i]*_xHit[i] + _yHit[i]*_yHit[i]);
    if ((Rz > Lower) && (Rz < Upper)) {
      double dZ = fabs(_zHit[i]-_zHit[i1]);
      if (dZ < dZmin) {
	dZmin = dZ;
	i3 = i;
      }
    }
  }

  // debug
  //std::cout << std::endl << Upper << "  " << Lower << "  " << dZmin << "  " << i3 << std::endl;




  double z1 = std::min(_zHit[i1],_zHit[i3]);
  double z3 = std::max(_zHit[i1],_zHit[i3]);


  int i2 = -1;
  double dRmin = 1.0e+20;
  double Rref = 0.5 * ( Rmax + Rmin );

  // 3d loop

  for (int i(0); i < _nHits; ++i) {
    if (_zHit[i] >= z1 && _zHit[i] <= z3) {
      double Rz = sqrt(_xHit[i]*_xHit[i] + _yHit[i]*_yHit[i]);
      double dRz = fabs(Rz - Rref);
      if (dRz < dRmin) {
	i2 = i;
	dRmin = dRz;
      }
    }
  }

  //int problematic = 0;

  if (i2 < 0 || i2 == i1 || i2 == i3) {
    //problematic = 1;
    // std::cout << "here we are " << std::endl;
    for (int i(0); i < _nHits; ++i) {
      if (i != i1 && i != i3) {
	i2 = i;
	if (_zHit[i2] < z1) {
	  int itemp = i1;
	  i1 = i2;
	  i2 = itemp;
	}
	else if (_zHit[i2] > z3) {
	  int itemp = i3;
	  i3 = i2;
	  i2 = itemp;
	}        
	break;
      }
    }      
    // std::cout << i1 << " " << i2 << " " << i3 << std::endl;
  }


  double x0  = 0.5*(_xHit[i2]+_xHit[i1]);
  double y0  = 0.5*(_yHit[i2]+_yHit[i1]);
  double x0p = 0.5*(_xHit[i3]+_xHit[i2]);
  double y0p = 0.5*(_yHit[i3]+_yHit[i2]);
  double ax  = _yHit[i2] - _yHit[i1];
  double ay  = _xHit[i1] - _xHit[i2];
  double axp = _yHit[i3] - _yHit[i2];
  double ayp = _xHit[i2] - _xHit[i3];
  double det = ax * ayp - axp * ay;
  double time;

  if (det == 0.) {
    time = 500.;
  }
  else {
    gsl_matrix* A = gsl_matrix_alloc(2,2);
    gsl_vector* B = gsl_vector_alloc(2);
    gsl_vector* T = gsl_vector_alloc(2);     
    gsl_matrix_set(A,0,0,ax);
    gsl_matrix_set(A,0,1,-axp);
    gsl_matrix_set(A,1,0,ay);
    gsl_matrix_set(A,1,1,-ayp);
    gsl_vector_set(B,0,x0p-x0);
    gsl_vector_set(B,1,y0p-y0);
    gsl_linalg_HH_solve(A,B,T);
    time = gsl_vector_get(T,0); 
    gsl_matrix_free(A);
    gsl_vector_free(B);
    gsl_vector_free(T);
  }

  double X0 = x0 + ax*time;
  double Y0 = y0 + ay*time;

  double dX = _xHit[i1] - X0;
  double dY = _yHit[i1] - Y0;

  double R0 = sqrt(dX*dX + dY*dY);

  /*
    if (problematic == 1) {
    std::cout << i1 << " " << i2 << " " << i3 << std::endl;
    std::cout << _xHit[i1] << " " << _yHit[i1] << " " << _zHit[i1] << std::endl;
    std::cout << _xHit[i2] << " " << _yHit[i2] << " " << _zHit[i2] << std::endl;
    std::cout << _xHit[i3] << " " << _yHit[i3] << " " << _zHit[i3] << std::endl;
    std::cout << "R0 = " << R0 << std::endl;
    }
  */

  double phi1 = (double)atan2(_yHit[i1]-Y0,_xHit[i1]-X0);
  double phi2 = (double)atan2(_yHit[i2]-Y0,_xHit[i2]-X0);
  double phi3 = (double)atan2(_yHit[i3]-Y0,_xHit[i3]-X0);

  // testing bz > 0 hypothesis

  if ( phi1 > phi2 ) 
    phi2 = phi2 + 2.0*M_PI;
  if ( phi1 > phi3 )
    phi3 = phi3 + 2.0*M_PI;
  if ( phi2 > phi3 )
    phi3 = phi3 + 2.0*M_PI;

  double bz_plus = (phi3 - phi1) / (_zHit[i3]-_zHit[i1]);
  double phi0_plus = phi1 - bz_plus * _zHit[i1];
  double dphi_plus = fabs( bz_plus * _zHit[i2] + phi0_plus - phi2 );

  // testing bz < 0 hypothesis

  phi1 = (double)atan2(_yHit[i1]-Y0,_xHit[i1]-X0);
  phi2 = (double)atan2(_yHit[i2]-Y0,_xHit[i2]-X0);
  phi3 = (double)atan2(_yHit[i3]-Y0,_xHit[i3]-X0);

  if ( phi1 < phi2 ) 
    phi2 = phi2 - 2.0*M_PI;
  if ( phi1 < phi3 )
    phi3 = phi3 - 2.0*M_PI;
  if ( phi2 < phi3 )
    phi3 = phi3 - 2.0*M_PI;

  double bz_minus = (phi3 - phi1) / (_zHit[i3]-_zHit[i1]);
  double phi0_minus = phi1 - bz_minus * _zHit[i1];
  double dphi_minus = fabs( bz_minus * _zHit[i2] + phi0_minus - phi2 );

  double bz;
  double phi0;

  if (dphi_plus < dphi_minus) {
    bz = bz_plus;
    phi0 = phi0_plus;
  }
  else {
    bz = bz_minus;
    phi0 = phi0_minus;

  }

  double par_init[5];

  if (parametrisation == 1) {
    par_init[0] = (double)X0;
    par_init[1] = (double)Y0;
    par_init[2] = (double)R0;
    par_init[3] = (double)bz;
    par_init[4] = (double)phi0;
  }
  else if (parametrisation == 2) {
    par_init[0] = (double)X0;
    par_init[1] = (double)Y0;
    par_init[2] = (double)(-phi0/bz);
    par_init[3] = (double)R0;
    par_init[4] = (double)(1/bz);
  }

  else if (parametrisation == 3) {  // parameter vector: (z0,Phi0,omega,d0,tanL)

    // debug
    // std::cout << std::setprecision(6) << "InitFit (X0,Y0,R0,bz,phi0) = " << "(" << X0 << "," << Y0 << "," << R0 << "," << bz << "," << phi0 << ")" << std::endl;

    
    // debug
    /*    
	  X0 = -1205.28;
	  Y0 = 175.317;
	  R0 = 1217.97;
	  bz = 0.00326074;
	  phi0 = -0.144444;
    */ 


    // debug
    // std::cout << std::setprecision(6) << "InitUsed (X0,Y0,R0,bz,phi0) = " << "(" << X0 << "," << Y0 << "," << R0 << "," << bz << "," << phi0 << ")" << std::endl;
        
    double omega = 1/R0*direction;
    double tanL  = (-1.0)*omega/bz;

    double Phi0  = (-1.0)*atan2(X0,Y0);
    
    if (direction == 1) {

      if (tanL >= 0.0) Phi0 += M_PI; // add pi (see LC-DET-2006-004) //  >= or > ?
      else Phi0 -= M_PI; // < or <= ?
      
    }

    //double d0 = R0 - X0/sin(Phi0);
    //double d0 = R0 + Y0/cos(Phi0);

    double d0 = 0.0;
    if (true /*direction != 1*/) d0 = R0 - sqrt(X0*X0 + Y0*Y0);
    // else d0 = R0 + sqrt(X0*X0 + Y0*Y0);



    // double d0 = R0 - ( (X0-Y0)/(sqrt(2.0)*cos(pi/4 - Phi0)) );    
    // double d0 = R0 - ((X0-Y0)/(sin(Phi0)+cos(Phi0)));

    // double Phi0 = asin(X0/(R0-d0));

    
    double z0 = (1/bz)*((-1.0)*phi0+Phi0+(omega*M_PI)/(2.0*fabs(omega)));


    // debug
    /*
      std::cout << std::setprecision(6) << "InitFitCalculated (d0,z0,phi0,omega,tanL) = " << "(" << d0 << "," << z0 << "," << Phi0 << "," << omega << "," << tanL << ")" 
      << "  " << "sign(omega) = " << direction << std::endl;
    */

    // debug        
    /*
      d0    = 0.00016512;
      z0    = 0.000853511;
      Phi0  = 1.11974;
      omega = -4.22171e-05;
      tanL  = -0.33436;
    */


    

    // debug
    // std::cout << std::setprecision(6) << "InitFitUsed (d0,z0,phi0,omega,tanL) = " << "(" << d0 << "," << z0 << "," << Phi0 << "," << omega << "," << tanL << ")" << std::endl;

    par_init[0] = z0;
    par_init[1] = Phi0;
    par_init[2] = omega;
    par_init[3] = d0;
    par_init[4] = tanL;
  }
  else return 1;


  // local variables
  int status = 0;
  int iter = 0;

  int npar = 5; // five parameters to fit
  int ndim = 0;
  if (parametrisation == 1) ndim = 2; // two dependent dimensions 
  else if (parametrisation == 2) ndim = 3; // three dependent dimensions
  else if (parametrisation == 3) ndim = 3; // three dependent dimensions

  else return 1;



  double chi2_nofit = 0.0;
  int iFirst = 1;
  for (int ipoint(0); ipoint < _nHits; ipoint++) {
    double distRPZ[2];
    double Dist = DistanceHelix(_xHit[ipoint],_yHit[ipoint],_zHit[ipoint],
				X0,Y0,R0,bz,phi0,distRPZ);
    double chi2rphi = distRPZ[0]/_exHit[ipoint];
    chi2rphi = chi2rphi*chi2rphi;
    double chi2z = distRPZ[1]/_ezHit[ipoint];
    chi2z = chi2z*chi2z;
    chi2_nofit = chi2_nofit + chi2rphi + chi2z;
    if (Dist > distmax || iFirst == 1) {
      distmax = Dist;
      iFirst = 0;
    }
  }      
  chi2_nofit = chi2_nofit/double(_nHits);

  if ( status_out == 1 ) {
    for (int i(0); i < 5; ++i) {
      parameter[i] = (double)par_init[i];
      dparameter[i] = 0.0;      
    }
    chi2 = chi2_nofit;
    return 0;
  }

  // converging criteria
  const double abs_error = 1e-4;
  const double rel_error = 1e-4;

  gsl_multifit_function_fdf fitfunct;

  const gsl_multifit_fdfsolver_type* T = gsl_multifit_fdfsolver_lmsder;

  gsl_multifit_fdfsolver* s = gsl_multifit_fdfsolver_alloc(T,ndim*_nHits,npar);

  gsl_matrix* covar = gsl_matrix_alloc(npar,npar);   // covariance matrix

  data d;
  d.n = _nHits;
  d.x = &_xHit[0];
  d.y = &_yHit[0];
  d.z = &_zHit[0];
  d.ex = &_exHit[0];
  d.ey = &_eyHit[0];
  d.ez = &_ezHit[0];


  if (parametrisation == 1) {
    fitfunct.f = &functParametrisation1;
    fitfunct.df = &dfunctParametrisation1;
    fitfunct.fdf = &fdfParametrisation1;
    fitfunct.n = ndim*_nHits;
    fitfunct.p = npar;
    fitfunct.params = &d;
  }
  else if (parametrisation == 2) {
    fitfunct.f = &functParametrisation2;
    fitfunct.df = &dfunctParametrisation2;
    fitfunct.fdf = &fdfParametrisation2;
    fitfunct.n = ndim*_nHits;
    fitfunct.p = npar;
    fitfunct.params = &d;
  }
  else if (parametrisation == 3) {
    fitfunct.f = &functParametrisation3;
    fitfunct.df = &dfunctParametrisation3;
    fitfunct.fdf = &fdfParametrisation3;
    fitfunct.n = ndim*_nHits;
    fitfunct.p = npar;
    fitfunct.params = &d;
  }
  else return 1;

  gsl_vector_view pinit = gsl_vector_view_array(par_init,npar);
  gsl_multifit_fdfsolver_set(s,&fitfunct,&pinit.vector);

  // perform fit
  do {
    iter++;
    status = gsl_multifit_fdfsolver_iterate(s);

    if (status) break;
    status = gsl_multifit_test_delta (s->dx, s->x,abs_error,rel_error);

  } while ( status==GSL_CONTINUE && iter < max_iter);

  //fg: jacobian has been dropped from gsl_multifit_fdfsolver in gsl 2:
  gsl_matrix * J = gsl_matrix_alloc(s->fdf->n, s->fdf->p);
  gsl_multifit_fdfsolver_jac( s, J);
  gsl_multifit_covar( J, rel_error, covar );
  //  gsl_multifit_covar (s->J, rel_error, covar);



  chi2 = 0.0;

  if (parametrisation == 1) {
    X0   = (double)gsl_vector_get(s->x,0);
    Y0   = (double)gsl_vector_get(s->x,1);
    R0   = (double)gsl_vector_get(s->x,2);
    bz   = (double)gsl_vector_get(s->x,3);
    phi0 = (double)gsl_vector_get(s->x,4);
  }
  else if (parametrisation == 2) {
    X0   = (double)gsl_vector_get(s->x,0);
    Y0   = (double)gsl_vector_get(s->x,1);
    R0   = (double)gsl_vector_get(s->x,3);
    bz   = (double)(1/gsl_vector_get(s->x,4));
    phi0 = (double)(-gsl_vector_get(s->x,2)/gsl_vector_get(s->x,4));
  }
  else if (parametrisation == 3) { // (parameter vector: (z0,phi0,omega,d0,tanL)

    double z0    = gsl_vector_get(s->x,0);
    double Phi0  = gsl_vector_get(s->x,1);
    double omega = gsl_vector_get(s->x,2);
    double d0    = gsl_vector_get(s->x,3);
    double tanL  = gsl_vector_get(s->x,4);

    X0   = (double)( ( (1/omega) - d0 )*sin(Phi0) );
    Y0   = (double)( (-1)*( (1/omega) - d0 )*cos(Phi0) );
    R0   = (double)( 1/fabs(omega) );
    bz   = (double)( (-1)*(omega/tanL) );
    phi0 = (double)( ( (z0*omega)/tanL ) + Phi0 + ( (omega*M_PI)/(2*fabs(omega)) ) );
  }
  else return 1;

  iFirst = 1;
  double ddmax = 0.0;
  for (int ipoint(0); ipoint < _nHits; ipoint++) {
    double distRPZ[2];
    double Dist = DistanceHelix(_xHit[ipoint],_yHit[ipoint],_zHit[ipoint],
				X0,Y0,R0,bz,phi0,distRPZ);
    double chi2rphi = distRPZ[0]/_exHit[ipoint];
    chi2rphi = chi2rphi*chi2rphi;
    double chi2z = distRPZ[1]/_ezHit[ipoint];
    chi2z = chi2z*chi2z;
    chi2 = chi2 + chi2rphi + chi2z;
    if (Dist > ddmax || iFirst == 1) {
      iFirst = 0;
      ddmax = Dist;
    }
  }
      

  chi2 = chi2/double(_nHits);
  if (chi2 < chi2_nofit) {
    for (int i = 0; i < npar; i++) {
      parameter[i]  = gsl_vector_get(s->x,i);
      dparameter[i] = sqrt(gsl_matrix_get(covar,i,i));
    }    
    distmax = ddmax;
  }
  else {
    chi2 = chi2_nofit;
    for (int i = 0; i < npar; i++) {
      parameter[i] = (double)par_init[i];
      dparameter[i] = 0.0;
    }
  }

  //  if (problematic == 1)
  //    std::cout << "chi2 = " << chi2 << std::endl;

  gsl_multifit_fdfsolver_free(s);
  gsl_matrix_free(covar);
  return 0; 

}



// ##########################################
// #####                                #####
// #####        private methods         #####
// #####                                #####
// ##########################################

//=============================================================================

void ClusterShapes::findElipsoid() {

  /**   Elipsoid parameter calculations see cluster_proper.f  */
  float cx,cy,cz ;
  float dx,dy,dz ;
  float r_hit_max, d_begn, d_last, r_max, proj;
  if (_ifNotInertia == 1) findInertia() ;
  //   Normalize the eigen values of inertia tensor
  float wr1 = sqrt(_ValAnalogInertia[0]/_totAmpl);
  float wr2 = sqrt(_ValAnalogInertia[1]/_totAmpl);
  float wr3 = sqrt(_ValAnalogInertia[2]/_totAmpl);
  _r1 = sqrt(wr2*wr3);                // spatial axis length -- the largest
  _r2 = sqrt(wr1*wr3);                // spatial axis length -- less
  _r3 = sqrt(wr1*wr2);                // spatial axis length -- even more less
  _vol = 4.*M_PI*_r1*_r2*_r3/3.;      // ellipsoid volume
  _r_ave = pow(_vol,1/3);             // average radius  (quibc root)
  _density = _totAmpl/_vol;           // density
  //    _eccentricity = _r_ave/_r1;   // Cluster Eccentricity
  _eccentricity =_analogWidth/_r1;   // Cluster Eccentricity

  // Find Minumal and Maximal Lenght for Principal axis
  r_hit_max = -100000.;
  d_begn    =  100000.;
  d_last    = -100000.;
  cx = _VecAnalogInertia[0] ;
  cy = _VecAnalogInertia[1] ;
  cz = _VecAnalogInertia[2] ;
  for (int i(0); i < _nHits; ++i) {
    dx = _xHit[i] - _xgr;
    dy = _yHit[i] - _ygr;
    dz = _zHit[i] - _zgr;
    r_max = sqrt(dx*dx + dy*dy + dz*dz);;
    if(r_max > r_hit_max) r_hit_max = r_max;
    proj = dx*cx + dy*cy + dz*cz;
    if(proj < d_begn)
      d_begn = proj;
    //            lad_begn = ladc(L)
    if(proj > d_last)
      d_last = proj;
    //            lad_last = ladc(L)
  }
  //        if (r_hit_max > 0.0)
  //	  _r1 = 1.05*r_hit_max; // + 5% of length
  _r1_forw = fabs(d_last);
  _r1_back = fabs(d_begn);
}

//=============================================================================

void ClusterShapes::findGravity() {

  _totAmpl = 0. ;
  for (int i(0); i < 3; ++i) {
    _analogGravity[i] = 0.0 ;
  }
  for (int i(0); i < _nHits; ++i) {
    _totAmpl+=_aHit[i] ;
    _analogGravity[0]+=_aHit[i]*_xHit[i] ;
    _analogGravity[1]+=_aHit[i]*_yHit[i] ;
    _analogGravity[2]+=_aHit[i]*_zHit[i] ;
  }
  for (int i(0); i < 3; ++i) {
    _analogGravity[i]/=_totAmpl ;
  }
  _xgr = _analogGravity[0];
  _ygr = _analogGravity[1];
  _zgr = _analogGravity[2];
  _ifNotGravity = 0;
}

//=============================================================================

void ClusterShapes::findInertia() {

  double aIne[3][3];
  //  float radius1;
  float radius2 = 0.0;

  findGravity();

  for (int i(0); i < 3; ++i) {
    for (int j(0); j < 3; ++j) {
      aIne[i][j] = 0.0;
    }
  }

  for (int i(0); i < _nHits; ++i) {
    float dX = _xHit[i] - _analogGravity[0];
    float dY = _yHit[i] - _analogGravity[1];
    float dZ = _zHit[i] - _analogGravity[2];
    aIne[0][0] += _aHit[i]*(dY*dY+dZ*dZ);
    aIne[1][1] += _aHit[i]*(dX*dX+dZ*dZ);
    aIne[2][2] += _aHit[i]*(dX*dX+dY*dY);
    aIne[0][1] -= _aHit[i]*dX*dY;
    aIne[0][2] -= _aHit[i]*dX*dZ;
    aIne[1][2] -= _aHit[i]*dY*dZ;
  }

  for (int i(0); i < 2; ++i) {
    for (int j = i+1; j < 3; ++j) {
      aIne[j][i] = aIne[i][j];
    }
  }
  //****************************************
  // analog Inertia
  //****************************************

  gsl_matrix_view aMatrix = gsl_matrix_view_array((double*)aIne,3,3);
  gsl_vector* aVector = gsl_vector_alloc(3);
  gsl_matrix* aEigenVec = gsl_matrix_alloc(3,3);
  gsl_eigen_symmv_workspace* wa = gsl_eigen_symmv_alloc(3);
  gsl_eigen_symmv(&aMatrix.matrix,aVector,aEigenVec,wa);
  gsl_eigen_symmv_free(wa);
  gsl_eigen_symmv_sort(aVector,aEigenVec,GSL_EIGEN_SORT_ABS_ASC);

  for (int i(0); i < 3; i++) {
    _ValAnalogInertia[i] = gsl_vector_get(aVector,i);
    for (int j(0); j < 3; j++) {
      _VecAnalogInertia[i+3*j] = gsl_matrix_get(aEigenVec,i,j);
    }
  }

  // Main principal points away from IP

  _radius = 0.;
  radius2 = 0.;

  for (int i(0); i < 3; ++i) {
    _radius += _analogGravity[i]*_analogGravity[i];
    radius2 += (_analogGravity[i]+_VecAnalogInertia[i])*(_analogGravity[i]+_VecAnalogInertia[i]);
  }

  if ( radius2 < _radius) {
    for (int i(0); i < 3; ++i)
      _VecAnalogInertia[i] = - _VecAnalogInertia[i];
  }

  _radius = sqrt(_radius);
  _ifNotInertia = 0;

  // The final job
  findWidth();
  findElipsoid();

  gsl_vector_free(aVector);
  gsl_matrix_free(aEigenVec);

}

//=============================================================================

void ClusterShapes::findWidth() {

  float dist = 0.0;
  if (_ifNotInertia == 1)  findInertia() ;
  _analogWidth  = 0.0 ;
  for (int i(0); i < _nHits; ++i) {
    dist = findDistance(i) ;
    _analogWidth+=_aHit[i]*dist*dist ;
  }
  _analogWidth  = sqrt(_analogWidth / _totAmpl) ;
  _ifNotWidth = 0 ;
}

//=============================================================================

float ClusterShapes::findDistance(int i) {

  float cx = 0.0;
  float cy = 0.0;
  float cz = 0.0;
  float dx = 0.0;
  float dy = 0.0;
  float dz = 0.0;
  cx = _VecAnalogInertia[0] ;
  cy = _VecAnalogInertia[1] ;
  cz = _VecAnalogInertia[2] ;
  dx = _analogGravity[0] - _xHit[i] ;
  dy = _analogGravity[1] - _yHit[i] ;
  dz = _analogGravity[2] - _zHit[i] ;
  float tx = cy*dz - cz*dy ;
  float ty = cz*dx - cx*dz ;
  float tz = cx*dy - cy*dx ;
  float tt = sqrt(tx*tx+ty*ty+tz*tz) ;
  float ti = sqrt(cx*cx+cy*cy+cz*cz) ;
  float f = tt / ti ;
  return f ;
}
/**
   Function sdist(xp,yp,zp,cx,cy,cz,xv,yv,zv)
   c----------------------------------------------------------------------
   c        Distance from line to point
   c       xp, yp, zp -- point is at the line
   c       xv, yv, zv -- point is out of line
   ********************************************************************
   *     Last update     V.L.Morgunov     08-Apr-2002                 *
   ********************************************************************
   real xp,yp,zp,cx,cy,cz,xv,yv,zv,t1,t2,t3,tt,sdist

   t1 = cy*(zp-zv)-cz*(yp-yv)
   t2 = cz*(xp-xv)-cx*(zp-zv)
   t3 = cx*(yp-yv)-cy*(xp-xv)
   tt = sqrt(cx**2+cy**2+cz**2)
   sdist = sqrt(t1**2+t2**2+t3**2)/tt

   return
   end
*/

//=============================================================================

float ClusterShapes::vecProduct(float * x1, float * x2) {

  float x1abs(0.);
  float x2abs(0.);
  float prod(0.);

  for (int i(0); i < 3; ++i) {
    x1abs += x1[i]*x1[i];
    x2abs += x2[i]*x2[i];
    prod  += x1[i]*x2[i];
  }


  x1abs = sqrt(x1abs);
  x2abs = sqrt(x2abs);

  if (x1abs > 0.0 && x2abs > 0.0) {
    prod = prod/(x1abs*x2abs);
  }
  else {
    prod = 0.;
  }

  return prod;

}

//=============================================================================

float ClusterShapes::vecProject(float * x, float * axis) {
  float axisabs(0.);
  float prod(0.);
  for (int i(0); i < 3; ++i) {
    axisabs += axis[i]*axis[i];
    prod  += x[i]*axis[i];
  }
  axisabs = sqrt(axisabs);
  if (axisabs > 0.0 ) {
    prod = prod/axisabs;
  }
  else {
    prod = 0.;
  }
  return prod;
}

//=============================================================================

double ClusterShapes::DistanceHelix(double x, double y, double z, double X0, double Y0,
				    double R0, double bz, double phi0, double * distRPZ) {

  double phi  = atan2(y-Y0,x-X0);
  double R    = sqrt( (y-Y0)*(y-Y0) + (x-X0)*(x-X0) );
  double dXY2 = (R-R0)*(R-R0);
  double _const_2pi = 2.0*M_PI;
  double xN = (bz*z + phi0 - phi)/_const_2pi;

  int n1 = 0;
  int n2 = 0;
  int nSpirals = 0;

  if (xN > 0) {
    n1 = (int)xN;
    n2 = n1 + 1;
  }
  else {
    n1 = (int)xN - 1;
    n2 = n1 + 1;
  }

  if (fabs(n1-xN) < fabs(n2-xN)) {
    nSpirals = n1;
  }
  else {
    nSpirals = n2;
  }

  double dZ = (phi + _const_2pi*nSpirals - phi0)/bz - z;
  double dZ2 = dZ*dZ;

  distRPZ[0] = sqrt(dXY2);
  distRPZ[1] = sqrt(dZ2);

  return sqrt(dXY2 + dZ2);

}

//=============================================================================

int ClusterShapes::transformToEigensystem(float* xStart, int& index_xStart, float* X0, float* Rm) {
  
  if (_ifNotInertia == 1) findInertia();

  float MainAxis[3];
  float MainCentre[3];


  MainAxis[0] = _VecAnalogInertia[0];
  MainAxis[1] = _VecAnalogInertia[1];
  MainAxis[2] = _VecAnalogInertia[2];
  MainCentre[0] = _analogGravity[0];
  MainCentre[1] = _analogGravity[1];
  MainCentre[2] = _analogGravity[2];

  //analoginertia looks something wrong!
  //change it to the direction of CoG
  //float ll=sqrt(MainCentre[0]*MainCentre[0]+MainCentre[1]*MainCentre[1]+MainCentre[2]*MainCentre[2]);
  //MainAxis[0]=MainCentre[0]/ll;
  //MainAxis[1]=MainCentre[1]/ll;
  //MainAxis[2]=MainCentre[2]/ll;

  int ifirst = 0;
  float xx[3];
  float prodmin = 1.0e+100;
  int index = 0;

  for (int i(0); i < _nHits; ++i) {
    xx[0] = _xHit[i] - MainCentre[0];
    xx[1] = _yHit[i] - MainCentre[1];
    xx[2] = _zHit[i] - MainCentre[2];
    float prod = vecProject(xx,MainAxis);
    if (ifirst == 0 || prod < prodmin) {
      ifirst = 1;
      prodmin = prod;
      index = i;
    }
  }
  xStart[0] = MainCentre[0] + prodmin*MainAxis[0];
  xStart[1] = MainCentre[1] + prodmin*MainAxis[1];
  xStart[2] = MainCentre[2] + prodmin*MainAxis[2];
  index_xStart = index;
  
  float l=sqrt(MainAxis[0]*MainAxis[0]+MainAxis[1]*MainAxis[1]+MainAxis[2]*MainAxis[2]);
  //float ll=sqrt(xStart[0]*xStart[0]+xStart[1]*xStart[1]+xStart[2]*xStart[2]);
  //std::cout << "xstart: " << index_xStart << " " << prodmin << " " << xStart[0] << std::endl;

  //calculate the surface of hcal
  float ecalrad=2.058000000e+03;   //in mm 
  float plugz=2.650000000e+03;   //in mm 
  
  float tmpcos=MainAxis[2]/l;
  float tmpsin=sqrt(MainAxis[0]*MainAxis[0]+MainAxis[1]*MainAxis[1])/l;
  float detend=0.0;
  if(fabs(xStart[2])<2.450000000e+03){   //if in the barrel
    detend=(ecalrad-sqrt(xStart[0]*xStart[0]+xStart[1]*xStart[1]))/tmpsin;
  }else{  //if in plug
    detend=(plugz-fabs(xStart[2]))/fabs(tmpcos);
  }
  if(detend<0.0) detend=0.0;

  // std::cout << "check length: " << detend << " "
  //           << xStart[0]/ll << " " << xStart[1]/ll << " " << xStart[2]/ll << " " 
  //           << MainAxis[0]/l << " " << MainAxis[1]/l << " " << MainAxis[1]/l << std::endl;
  
  for (int i(0); i < _nHits; ++i) {
    xx[0] = _xHit[i] - xStart[0];
    xx[1] = _yHit[i] - xStart[1];
    xx[2] = _zHit[i] - xStart[2];
    float xx2(0.);
    for (int j(0); j < 3; ++j) xx2 += xx[j]*xx[j];
    
    _xl[i] = 0.001 + vecProject(xx,MainAxis);
    _xt[i] = sqrt(std::max(0.0,xx2 + 0.01 - _xl[i]*_xl[i]));
    //    std::cout << i << " " << _xl[i] << " " << _xt[i] << " " << _aHit[i] << " "
    //              << std::endl;
  }
  
  //first, check ecal and solve wrong behaviour
  //std::cout << "param: " << X0[0] << " " << X0[1] << " " << Rm[0] << " " << Rm[1] << std::endl;
  for (int i = 0; i < _nHits; ++i) { 
    //if(_types[i]==0 || _types[i]==3){   //if the hit is in Ecal
    _t[i] = _xl[i]/X0[0];
    _s[i] = _xt[i]/Rm[0];
    if(detend<_xl[i]){
      //std::cout << "something wrong!! " << detend << " " << _xl[i] << " " << _t[i] << std::endl;
      //detend = _xl[i]; //to avoid wrong behaviour
    }
    //}
  }
  
  //second, check hcal
  /*for (int i = 0; i < _nHits; ++i) { 
    if(_types[i]==1 || _types[i]==4){   //if the hit is in Hcal
    if(_xl[i]>detend){
    _t[i] = detend/X0[0]+(_xl[i]-detend)/X0[1];
    _s[i] = _xt[i]/Rm[1];
    }else{
    _t[i] = _xl[i]/X0[0];
    _s[i] = _xt[i]/Rm[0];
    }   
    }
    // std::cout << _types[i] << " " << _xl[i] << " "
    //           << _xl[i]+sqrt(xStart[0]*xStart[0]+xStart[1]*xStart[1]+xStart[2]*xStart[2])  << " "
    //           << _t[i] << " " << _s[i] << std::endl;
    }*/
  
  _ifNotEigensystem = 0;
  
  return 0; // no error messages at the moment

}

//=============================================================================

float ClusterShapes::calculateChi2Fit3DProfileSimple(float a, float b, float c,
						     float d) {

  // ClusterShapes::transformToEigensystem needs to be executed before

  float chi2 = 0.0;
  float Ampl = 0.0; 

  for (int i(0); i < _nHits; ++i) {
    // old definition of Ampl and chi2
    //Ampl = a*(float)pow(_xl[i],b)*exp(-c*_xl[i]-d*_xt[i]); 
    //chi2 += ((Ampl - _aHit[i])*(Ampl - _aHit[i]))/(_aHit[i]*_aHit[i]);

    Ampl = a*(float)pow(_xl[i],b)*exp(-c*_xl[i]-d*_xt[i]);
    chi2 += (log(_aHit[i]) - log(Ampl))*(log(_aHit[i]) - log(Ampl));
  }

  chi2 = chi2/std::max((float)1.0,(float)(_nHits - 4));

  return chi2;

}

//=============================================================================

float ClusterShapes::calculateChi2Fit3DProfileAdvanced(float E0, float a, float b,
						       float d, float t0) {

  // ClusterShapes::transformToEigensystem needs to be executed before

  float chi2  = 0.0;
  float Ampl  = 0.0;
  float shift = 0.0;

  for (int i(0); i < _nHits; ++i) {

    shift = _t[i]-t0;

    if (shift <= 0) Ampl = 0.0;
    else {
      Ampl = E0 * b * invG(a) * pow(b*(shift),a-1) 
	* exp(-b*(shift)) * exp(-d*_s[i]);
    }

    // debug 
    /*
      std::cout << "OUT : " << Ampl << "  " << E0  << "  " << a  << "  " << b  << "  " 
      << d  << "  " << t0  << "  " << _t[i] << "  " << invG(a)  << "  "
      << pow(b*(_t[i]-t0),a-1)  << "  " << exp(-b*(_t[i]-t0))  << "  " 
      << exp(-d*_s[i]) 
      << std::endl;
    */

    chi2 += ((Ampl - _aHit[i])*(Ampl - _aHit[i]))/(_aHit[i]*_aHit[i]);
    //chi2 += (log(_aHit[i]) - log(Ampl))*(log(_aHit[i]) - log(Ampl));
    //chi2 += log((Ampl - _aHit[i])*(Ampl - _aHit[i]))/log(_aHit[i]*_aHit[i]);

  }
  chi2 = chi2/std::max((float)1.0,(float)(_nHits - 4));

  return chi2;

}

//=============================================================================

int ClusterShapes::fit3DProfileSimple(float& chi2, float& a, float& b, float& c,
				      float& d) {

  // Fit function: _aHit[i] = a * pow(_xl[i],b) * exp(-c*_xl[i]) * exp(-d*_xt[i])

  float Slnxl(0.);
  float Sxl(0.);
  float Sxt(0.);
  float Sln2xl(0.);
  float Sxllnxl(0.);
  float Sxtlnxl(0.);
  float Sxlxl(0.);
  float Sxlxt(0.);
  float Sxtxt(0.);
  float SlnA(0.);
  float SlnAlnxl(0.);
  float SlnAxl(0.);
  float SlnAxt(0.);

  // for a quadratic matrix
  for (int i = 0; i < _nHits; i++) {
    Slnxl += log(_xl[i]);
    Sxl += _xl[i];
    Sxt += _xt[i];
    Sln2xl += log(_xl[i])*log(_xl[i]);
    Sxllnxl += _xl[i]*log(_xl[i]);
    Sxtlnxl += _xt[i]*log(_xl[i]);
    Sxlxl += _xl[i]*_xl[i];
    Sxlxt += _xl[i]*_xt[i];
    Sxtxt += _xt[i]*_xt[i];
    SlnA += log(_aHit[i]);
    SlnAlnxl += log(_aHit[i])*log(_xl[i]);
    SlnAxl += log(_aHit[i])*_xl[i];
    SlnAxt += log(_aHit[i])*_xt[i]; 
  }
  // create system of linear equations, written as Ae = z

  gsl_matrix* A = gsl_matrix_alloc(4,4);
  gsl_vector* z = gsl_vector_alloc(4);
  gsl_vector* e = gsl_vector_alloc(4);
  
  // initialise matrix and vectors
  
  gsl_matrix_set(A,0,0,_nHits);
  gsl_matrix_set(A,0,1,Slnxl);
  gsl_matrix_set(A,0,2,-Sxl);
  gsl_matrix_set(A,0,3,-Sxt);
  
  gsl_matrix_set(A,1,0,Slnxl);
  gsl_matrix_set(A,1,1,Sln2xl);
  gsl_matrix_set(A,1,2,-Sxllnxl);
  gsl_matrix_set(A,1,3,-Sxtlnxl);
  
  gsl_matrix_set(A,2,0,-Sxl);
  gsl_matrix_set(A,2,1,-Sxllnxl);
  gsl_matrix_set(A,2,2,Sxlxl);
  gsl_matrix_set(A,2,3,Sxlxt);

  gsl_matrix_set(A,3,0,-Sxt);
  gsl_matrix_set(A,3,1,-Sxtlnxl);
  gsl_matrix_set(A,3,2,Sxlxt);
  gsl_matrix_set(A,3,3,Sxtxt);

  gsl_vector_set(z,0,SlnA);
  gsl_vector_set(z,1,SlnAlnxl);
  gsl_vector_set(z,2,-SlnAxl);
  gsl_vector_set(z,3,-SlnAxt);

  gsl_linalg_HH_solve(A,z,e);

  a = exp(gsl_vector_get(e,0));
  b = gsl_vector_get(e,1);
  c = gsl_vector_get(e,2);
  d = gsl_vector_get(e,3);

  chi2 = calculateChi2Fit3DProfileSimple(a,b,c,d);
  
  gsl_matrix_free(A);
  gsl_vector_free(z);
  gsl_vector_free(e);

  int result = 0;  // no error handling at the moment
  return result;

}

//=============================================================================

int ClusterShapes::fit3DProfileAdvanced(float& chi2, double* par_init, double* par,
					int npar, float* t, float* s, float* E, 
					float E0) {

  // local variables
  int status = 0;
  int iter = 0;
  int max_iter = 1000;


  // converging criteria
  const double abs_error = 0.0;
  const double rel_error = 1e-1;



  gsl_multifit_function_fdf fitfunct;

  const gsl_multifit_fdfsolver_type* T = gsl_multifit_fdfsolver_lmsder;

  //std::cout << T << " " << _nHits << " " << npar << std::endl;

  gsl_multifit_fdfsolver* Solver = gsl_multifit_fdfsolver_alloc(T,_nHits,npar);

  gsl_matrix* covar = gsl_matrix_alloc(npar,npar);   // covariance matrix

  data DataSet;
  DataSet.n = _nHits;
  DataSet.x = &t[0];
  DataSet.y = &s[0];
  DataSet.z = &E[0];  // _aHit[0]; // ???? normalise per volume ????


  fitfunct.f = &ShapeFitFunct;
  fitfunct.df = &dShapeFitFunct;
  fitfunct.fdf = &fdfShapeFitFunct;
  fitfunct.n = _nHits;
  fitfunct.p = npar;
  fitfunct.params = &DataSet;

  gsl_vector_view pinit = gsl_vector_view_array(par_init,npar);

  //std::cout << Solver << " " << &fitfunct << " " << &pinit.vector << std::endl;
  gsl_multifit_fdfsolver_set(Solver,&fitfunct,&pinit.vector);

  gsl_set_error_handler_off();

  // perform fit
  do {
    iter++;
    //std::cout << "Multidimensional Fit Iteration started ... ... ";
    status = gsl_multifit_fdfsolver_iterate(Solver);
    //std::cout << "--- DONE ---" << std::endl;

    if (status) break;
    status = gsl_multifit_test_delta (Solver->dx,Solver->x,abs_error,rel_error);

    // debug
    /*
    //    E0  = (float)gsl_vector_get(Solver->x,0);
    par[0] = (float)gsl_vector_get(Solver->x,0);
    par[1] = (float)gsl_vector_get(Solver->x,1);
    par[2] = (float)gsl_vector_get(Solver->x,2);
    par[3] = (float)gsl_vector_get(Solver->x,3);

    std::cout << "Status of multidimensional fit : " << status << "  "
    << "Iterations : " << iter << std::endl;
    std::cout << "E0 : " <<  "FIXED" << "\t" << "A : " << par[0] << "\t" << "B : " 
    <<  par[1] << "\t" << "D : " << par[2] << "\t" << "t0 : "  << par[3] 
    << std::endl << std::endl;
    */

  } while ( status==GSL_CONTINUE && iter < max_iter);

  //fg: jacobian has been dropped from gsl_multifit_fdfsolver in gsl 2:
  gsl_matrix * J = gsl_matrix_alloc(Solver->fdf->n, Solver->fdf->p);
  gsl_multifit_fdfsolver_jac(Solver, J);
  gsl_multifit_covar( J, rel_error, covar );
  //  gsl_multifit_covar (Solver->J,rel_error,covar);

  //  E0  = (float)gsl_vector_get(Solver->x,0);
  par[0] = (float)gsl_vector_get(Solver->x,0); // A
  par[1] = (float)gsl_vector_get(Solver->x,1); // B
  par[2] = (float)gsl_vector_get(Solver->x,2); // D
  par[3] = (float)gsl_vector_get(Solver->x,3); // t0

  gsl_multifit_fdfsolver_free(Solver);
  gsl_matrix_free(covar);
  
  chi2 = calculateChi2Fit3DProfileAdvanced(E0,par[0],par[1],par[2],par[3]);
  if (status) chi2 = -1.0;

  int result = 0;  // no error handling at the moment
  return result;

}

//=============================================================================

float ClusterShapes::getEmax(float* xStart, int& index_xStart, float* X0, float* Rm){

  if (_ifNotEigensystem == 1){
    transformToEigensystem(xStart,index_xStart,X0,Rm);
  }

  float E_max=0.0;
  //int i_max=0;
  for (int i = 0; i < _nHits; ++i) {
    if (E_max < _aHit[i]) {
      E_max = _aHit[i]; 
      //i_max = i;
    }
  }
  
  return E_max;
}

float ClusterShapes::getsmax(float* xStart, int& index_xStart, float* X0, float* Rm){

  if (_ifNotEigensystem == 1){
    transformToEigensystem(xStart,index_xStart,X0,Rm);
  }

  float E_max=0.0,xl_max=0.0;
  float xl_start=1.0e+50;
  //int i_max=0;
  for (int i = 0; i < _nHits; ++i) {
    //check the position of maximum energy deposit
    if (E_max < _aHit[i]) {
      E_max = _aHit[i]; 
      xl_max = _xl[i];
      //i_max = i;
    }
    //check the position of shower start
    if (xl_start > _xl[i]) {
      xl_start = _xl[i];
    }
  }
  
  return fabs(xl_max-xl_start);
}

float ClusterShapes::getxl20(float* xStart, int& index_xStart, float* X0, float* Rm){
  
  if (_ifNotEigensystem == 1){
    transformToEigensystem(xStart,index_xStart,X0,Rm);
  }
  
  float *xl_res = new float[_nHits];
  float *E_res = new float[_nHits];
  float E_tot=0.0;
  for (int i = 0; i < _nHits; ++i) {
    E_res[i]=_aHit[i];
    xl_res[i]=_xl[i];

    E_tot+=_aHit[i];
  }
  
  //sort deposite energy to xs ascending order
  float tmpe=0.0, tmpxl=0.0;
  for (int i = 0; i < _nHits; ++i) {
    for (int j = i + 1; j < _nHits; ++j) {
      if(xl_res[i]>xl_res[j]){
	tmpe=E_res[i];
	E_res[i]=E_res[j];
	E_res[j]=tmpe;

	tmpxl=xl_res[i];
	xl_res[i]=xl_res[j];
	xl_res[j]=tmpxl;
      }
    }
  }

  // std::cout << "check xt: " << xt_res[0] << " " << xt_res[1] << " " << xt_res[2]
  //           <<std::endl;
  float E20=0.0,okxl=0.0;
  int k=0;
  while(E20/E_tot<0.2){
    E20+=E_res[k];
    k++;
  }
  
  //final hit is located in outer radius
  okxl=xl_res[k-2];

  delete[] xl_res;
  delete[] E_res;

  return okxl;
}

float ClusterShapes::getxt90(float* xStart, int& index_xStart, float* X0, float* Rm){
  
  if (_ifNotEigensystem == 1){
    transformToEigensystem(xStart,index_xStart,X0,Rm);
  }
  
  float *xt_res = new float[_nHits];
  float *E_res = new float[_nHits];
  float E_tot=0.0;
  for (int i = 0; i < _nHits; ++i) {
    E_res[i]=_aHit[i];
    xt_res[i]=_xt[i];

    E_tot+=_aHit[i];
  }
  
  //sort deposite energy to xt ascending order
  float tmpe=0.0, tmpxt=0.0;
  for (int i = 0; i < _nHits; ++i) {
    for (int j = i + 1; j < _nHits; ++j) {
      if(xt_res[i]>xt_res[j]){
	tmpe=E_res[i];
	E_res[i]=E_res[j];
	E_res[j]=tmpe;

	tmpxt=xt_res[i];
	xt_res[i]=xt_res[j];
	xt_res[j]=tmpxt;
      }
    }
  }

  // std::cout << "check xt: " << xt_res[0] << " " << xt_res[1] << " " << xt_res[2]
  //           <<std::endl;
  float E90=0.0,okxt=0.0;
  int k=0;
  while(E90/E_tot<0.9){
    E90+=E_res[k];
    k++;
  }
  
  //final hit is located in outer radius
  okxt=xt_res[k-2];

  delete[] xt_res;
  delete[] E_res;

  return okxt;
}

//for test
void ClusterShapes::gethits(float* xStart, int& index_xStart, float* X0, float* Rm, float *okxl, float *okxt, float *oke){
  
  if (_ifNotEigensystem == 1){
    transformToEigensystem(xStart,index_xStart,X0,Rm);
  }
  
  for (int i = 0; i < _nHits; ++i) {
    okxl[i]=_xl[i];
    okxt[i]=_xt[i];
    oke[i]=_aHit[i];
  }
  
  return;
}


float ClusterShapes::getRhitMean(float* xStart, int& index_xStart, float* X0, float* Rm){

  if (_ifNotEigensystem == 1){
    transformToEigensystem(xStart,index_xStart,X0,Rm);
  }

  float MainCentre[3];

  MainCentre[0] = _analogGravity[0];
  MainCentre[1] = _analogGravity[1];
  MainCentre[2] = _analogGravity[2];

  // float Rhit[_nHits]={0};
  float Rhit=0;

  float Rhitsum=0;
  float Rhitmean=0;

  for (int i = 0; i < _nHits; ++i) {
    Rhit = sqrt(pow((_xHit[i]-MainCentre[0]),2) + pow((_yHit[i]-MainCentre[1]),2));
    Rhitsum += Rhit;
  }
  // cogx-- xx[0]  
  
  Rhitmean = Rhitsum/_nHits;

  return Rhitmean;
}

float ClusterShapes::getRhitRMS(float* xStart, int& index_xStart, float* X0, float* Rm){

  if (_ifNotEigensystem == 1){
    transformToEigensystem(xStart,index_xStart,X0,Rm);
  }

  float MainCentre[3];

  MainCentre[0] = _analogGravity[0];
  MainCentre[1] = _analogGravity[1];
  MainCentre[2] = _analogGravity[2];

  // float Rhit[_nHits]={0};
  float Rhit=0;

  float Rhit2sum=0;
  float Rhitrms=0;

  for (int i = 0; i < _nHits; ++i) {
    Rhit = sqrt(pow((_xHit[i]-MainCentre[0]),2) + pow((_yHit[i]-MainCentre[1]),2));
    Rhit2sum += pow(Rhit,2);
  }
  // cogx-- xx[0]  
  
  Rhitrms = sqrt(Rhit2sum/_nHits);

  return Rhitrms;
}


