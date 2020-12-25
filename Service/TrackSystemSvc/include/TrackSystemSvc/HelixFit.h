#ifndef __HELIXFIT__
#define __HELIXFIT__
//
//  HelixFit.h
//  MarlinTrk
//

/**
 //  Created by Steve Aplin on 9/16/11.
 //  DESY
 //
 //  C++ rewrite of the aleph Fortran routine TFITHL
 //
 //! Fast helix fit 
 //
 //
 //   Input:  NPT           Number of 3-D points to be fit
 //           xf       Array of X-values of points to be fit
 //           yf       Array of Y-values of points to be fit
 //           zf       Array of Z-values of points to be fit
 //           wf       Array of 1/(sig(rphi))**2 for each point
 //           wzf      Array of 1/(sig(z))**2 for each point
 //           iopt     < 3 : error matrix calculated 
 //                    = 3 : 3-dimensional iteration
 //
 //  OUTPUT: vv0    = Helix parameter in perigee form
 //          ee0    = INVERSE OF ERROR MATRIX IN TRIANG. FORM
 //          chi2ph = CHI SQUARED = SUM (PHI DEVIATIONS/ERRORS)**2
 //          CH2Z  = CHI SQUARED = SUM (Z DEVIATIONS/ERRORS)**2
 //  NOTE: DEGREES OF FREEDOM = 2*NPT-5
 //----------------------------------------------------------------
 //     BASED ON  SUBROUTINE CIRCLE
 //     REFERENCE:  COMPUTER PHYSICS COMMUNICATIONS VOL 33,P329
 //
 //   AUTHORS:  N. CHERNOV, G. OSOSKOV & M. POPPE
 //   Modified by:  Fred Weber, 8 Jun 1989
 //   Modified by:  M.Cattaneo, 27-Jan-1998
 //                 Protect against arg SIN > 1.0
 //
 //-----------------------------------------------------------------
 */

namespace MarlinTrk {
  
  class HelixFit {
    
    
  public:
    
    int fastHelixFit(int npt, double* xf, double* yf, float* rf, float* pf, double* wf, float* zf , float* wzf, int iopt,
                     float* vv0, float* ee0, float& ch2ph, float& ch2z);
    
  };
  
}

#endif
