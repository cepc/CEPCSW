#ifndef _Arbor_hh_
#define _Arbor_hh_

#include "ArborHit.h"
#include <string>
#include <iostream>
#include <TVector3.h>
#include "ArborTool.h"

void init();

void HitsCleaning( std::vector<ArborHit> inputHits );

void HitsClassification( linkcoll inputLinks );

void BuildInitLink(std::vector<float> Thresholds);

void LinkIteration(int time);

void BranchBuilding(float SeedThreshold);

branchcoll Arbor( std::vector<ArborHit>, std::vector<float> Thresholds );

/*
 int NLayer_A, NStave_A, SubD_A; 
 int NLayer_B, NStave_B, SubD_B;
 float MagA, MagB, Depth_A, Depth_B, ECCorr, DisAB; 
 int FlagTrkPS, FlagEH; 
 int FlagPSEE, FlagHH;       
 int FlagStaveSame;
 int FlagStaveDiff;
 TVector3 PosA, PosB, PosDiffAB, PosDiffBA, linkDir;
*/

#endif


