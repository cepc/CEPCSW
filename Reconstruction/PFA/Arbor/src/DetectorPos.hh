#ifndef DETECTORPOS_H_
#define DETECTORPOS_H_

#include "TVector3.h"
#include "TMatrixF.h"
#include <iostream>
#include <vector>
#include <string>

const float BField = 3.0;
const float DHCALBarrelRadius = 2058.0;         //Octo
const float DHCALEndCapInnerZ = 2650.0;
const float ECALEndCapInnerZ = 2450.0;
const float ECALHalfZ = 2350.0; // mm, Endcap Ente
const float ECALRadius = 1847.4; // mm... minimal part for the octagnle.
const float TPCRadius = 1808.0 ;

const float TPCOuterRadius = 1808.0; 
const float TPCInnerRadius = 325.0; 
const float LStar = 1500.0; 

const float DHCALCalibrationConstant = 0.12;

int BarrelFlag( TVector3 inputPos );

int DepthFlag( TVector3 inputPos );

TVector3 CalVertex( TVector3 Pos1, TVector3 Dir1, TVector3 Pos2, TVector3 Dir2 );

int TPCPosition( TVector3 inputPos );		//Used to tag MCParticle position, if generated inside TPC & Dead outside

float DisSeedSurface( TVector3 SeedPos );	//for a given position, calculate the distance to Calo surface ( ECAL )

float DisSeedSurfaceSimple( TVector3 SeedPos );


float DisTPCBoundary( TVector3 Pos );

float DistanceChargedParticleToCluster(TVector3 CPRefPos, TVector3 CPRefMom, TVector3 CluPosition);

#endif //
