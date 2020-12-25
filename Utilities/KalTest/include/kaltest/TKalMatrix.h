#ifndef TKALMATRIX_H
#define TKALMATRIX_H
//*************************************************************************
//* ===================
//*  TKalMatrix Class
//* ===================
//*
//* (Description)
//*    TKalMatrix is a wrapper of TMatrixD.
//* (Requires)
//* 	TMatrixD
//* (Provides)
//* 	class TKalMatrix
//* (Update Recored)
//*   2003/09/30  K.Fujii	Original version.
//*
//*************************************************************************

#include "TMatrixD.h"
#include "TVector3.h"
//_____________________________________________________________________
//  ------------------------------
//  Base Class for matrix used by Kalman filter
//  ------------------------------
//
class TKalMatrix : public TMatrixD { 
public:
   TKalMatrix(Int_t rowdim = 1, Int_t coldim = 1);

   TKalMatrix(const TKalMatrix &orig);
   TKalMatrix(const TMatrixD &orig);
                                                                                
   TKalMatrix(TMatrixD::EMatrixCreatorsOp1 op,
              const TKalMatrix &prototype);
   TKalMatrix(TMatrixD::EMatrixCreatorsOp1 op,
              const TMatrixD &prototype);
                                                                                
   TKalMatrix(const TKalMatrix &a,
              TMatrixD::EMatrixCreatorsOp2 op,
              const TKalMatrix &b) ;
   TKalMatrix(const TMatrixD &a,
              TMatrixD::EMatrixCreatorsOp2 op,
              const TMatrixD &b) ;

   TKalMatrix(const TVector3 &v);

   virtual ~TKalMatrix() {}

   virtual void      DebugPrint(Option_t *opt = "", Int_t nc = 5) const;

   static TKalMatrix ToKalMat  (const TVector3 &vec);
   static TVector3   ToThreeVec(const TMatrixD &mat);
   
private:
   
   ClassDef(TKalMatrix,1)      // Base class for Kalman matrix
};
#endif
