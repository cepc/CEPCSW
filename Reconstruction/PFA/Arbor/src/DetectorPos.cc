#include <TMath.h>
#include "DetectorPos.hh"

using namespace std;

const double pi = acos(-1.0); 

//Geometric Parameter - ... need to be changed for different detector models

int BarrelFlag( TVector3 inputPos)
{
	int isBarrel = 0;	

	if(fabs(inputPos[2]) < ECALHalfZ)
	{
		isBarrel = 1; 
	}

	return isBarrel; 
}

int DepthFlag( TVector3 inputPos )      //Used to calculate depth of given position...
{

	float ShiftDHCAL = 530; // 20 Layers Inside DHCAL

	float DHCALDeepInnerZ = DHCALEndCapInnerZ + ShiftDHCAL;
	float DHCALDeepBarrelRadius = DHCALBarrelRadius + ShiftDHCAL;

	float DHCALInnerOctRadius = DHCALBarrelRadius/cos(pi/4.0);
	float DHCALDeepOctRadius = DHCALBarrelRadius/cos(pi/4.0) + ShiftDHCAL;
	float ECALInnerOctRadius = ECALRadius/cos(pi/4.0);

	int FlagD(-1);

	if( fabs(inputPos[2]) > DHCALDeepInnerZ || fabs(inputPos[1]) > DHCALDeepBarrelRadius || fabs(inputPos[0]) > DHCALDeepBarrelRadius || fabs(inputPos[0] + inputPos[1]) > DHCALDeepOctRadius || fabs(inputPos[0] - inputPos[1]) > DHCALDeepOctRadius )
	{
		FlagD = 2;
	}
	else if( fabs(inputPos[2]) > DHCALEndCapInnerZ || fabs(inputPos[1]) > DHCALBarrelRadius || fabs(inputPos[0]) > DHCALBarrelRadius || fabs(inputPos[0] + inputPos[1]) > DHCALInnerOctRadius || fabs(inputPos[0] - inputPos[1]) > DHCALInnerOctRadius )
	{
		FlagD = 1;          // Position outsider than DHCAL Region
	}
	else if( fabs(inputPos[2]) > ECALEndCapInnerZ || fabs(inputPos[1]) > ECALRadius || fabs(inputPos[0]) > ECALRadius || fabs(inputPos[0] + inputPos[1]) > ECALInnerOctRadius || fabs(inputPos[0] - inputPos[1]) > ECALInnerOctRadius )
	{
		FlagD = 0;
	}
	else
	{
		FlagD = 10;         // Position inside Calo... Problematic for Seeds... But could be PreShower hits.
	}

	return FlagD;

}

TVector3 CalVertex( TVector3 Pos1, TVector3 Dir1, TVector3 Pos2, TVector3 Dir2 )
{
	TVector3 VertexPos;

	float tau1(0), tau2(0);

	TVector3 delP;
	delP = Pos1 - Pos2;

	double Normal(0);
	Normal = (Dir1.Dot(Dir2))*(Dir1.Dot(Dir2)) - Dir1.Mag()*Dir1.Mag()*Dir2.Mag()*Dir2.Mag();

	if(Normal != 0)
	{
		tau1 = (Dir2.Mag()*Dir2.Mag() * delP.Dot(Dir1) - Dir1.Dot(Dir2)*delP.Dot(Dir2))/Normal;
		tau2 = (Dir1.Dot(Dir2)*delP.Dot(Dir1) - Dir1.Mag()*Dir1.Mag() * delP.Dot(Dir2))/Normal;
	}

	VertexPos = 0.5*(Pos1 + Pos2 + tau1*Dir1 + tau2*Dir2);

	return VertexPos;
}

int TPCPosition( TVector3 inputPos )
{
	int flagPos(-1); // == 0 means inside TPC, == 1 means outside; 

	if( fabs(inputPos[2]) > ECALHalfZ || sqrt( inputPos[0]*inputPos[0] + inputPos[1]*inputPos[1] ) > TPCRadius ) flagPos = 1;
	else flagPos = 0;

	return flagPos;
}

float DisSeedSurface( TVector3 SeedPos )	//ECAL, HCAL, EndCapRing...
{

	float DisSS = 0;

	if( fabs(SeedPos[2]) > ECALHalfZ )         //EcalEndcap hit start from 2350 + 100 = 2450
	{

		if( SeedPos.Perp() > ECALRadius )
		{
			if( fabs(SeedPos[2])/SeedPos.Perp() > (ECALHalfZ + 103)/(ECALRadius + 100) )
			{
				DisSS = ( fabs(SeedPos[2]) - ECALHalfZ - 103 ) * SeedPos.Mag()/fabs(SeedPos[2]);
			}
			else
			{
				DisSS = (SeedPos.Perp() - ECALRadius - 100 )*SeedPos.Mag()/SeedPos.Perp();
			}
		}
		else
		{
			DisSS = fabs(SeedPos[2]) - ECALHalfZ - 103;
		}
	}
	else if( SeedPos.Perp() > ECALRadius + 400 )
	{
		DisSS = SeedPos.Perp() - ECALRadius - 100;
	}
	else if( (SeedPos.Phi() > 0 && int(SeedPos.Phi() * 4/pi + 0.5) % 2 == 0 ) || (SeedPos.Phi() < 0 && int(SeedPos.Phi() * 4/pi + 8.5) % 2 == 0 ))
	{
		DisSS = min( fabs(fabs(SeedPos[0]) - ECALRadius), fabs(fabs(SeedPos[1]) - ECALRadius ) );
	}
	else
	{
		DisSS = min( fabs(fabs(SeedPos[0] + SeedPos[1])/1.414214 -ECALRadius), fabs(fabs(SeedPos[0] - SeedPos[1])/1.414214 - ECALRadius) );
	}

	return DisSS;
}

float DisSeedSurfaceSimple( TVector3 SeedPos )        //ECAL, HCAL, EndCapRing...
{
	// TVector3 SeedPos = a_clu->getPosition();
	float DisSS = 0;
	float DisZ = abs(SeedPos.Z()) - ECALHalfZ;
	float DisR = SeedPos.Perp() - ECALRadius;
	if(DisR < 0 && DisZ > 0)
	{
		DisSS = DisZ; 
	}
	else if(DisZ < 0 && DisR > 0)
	{
		DisSS = DisR;
	}
	else if(DisR > 0 && DisZ > 0)
	{
		if( fabs(SeedPos.Z())/SeedPos.Perp() > ECALHalfZ/ECALRadius)
		{
			DisSS = ( fabs(SeedPos[2]) - ECALHalfZ) * SeedPos.Mag()/fabs(SeedPos[2]);
		}
		else
		{
			DisSS = (SeedPos.Perp() - ECALRadius)*SeedPos.Mag()/SeedPos.Perp();
		}
	}

	return DisSS;
}

/*
float DisSeedSurfaceClu( Cluster * a_clu )        //ECAL, HCAL, EndCapRing...
{
	TVector3 SeedPos = a_clu->getPosition();
	return DisSeedSurface(SeedPos);
}
*/

float DisTPCBoundary( TVector3 Pos )
{
	float DisZ = TMath::Min( fabs(ECALHalfZ-Pos.Z()),fabs(ECALHalfZ+Pos.Z()) );
	float DisR = ECALRadius - Pos.Perp();  
	float Dis = TMath::Min(DisZ, DisR);

	return Dis; 
}

float DistanceChargedParticleToCluster(TVector3 CPRefPos, TVector3 CPRefMom, TVector3 CluPosition)	//Extend to Track/MCP
{
	// Line extrapolation from RefPos with RefMom, calculate the minimal distance to Cluster

	float DisCPClu = 0; 
	TVector3 Diff_Clu_CPRef, NormCPRefMom; 

	Diff_Clu_CPRef = CluPosition - CPRefPos; 
	NormCPRefMom = 1.0/CPRefMom.Mag()*CPRefMom;
	float ProDis = Diff_Clu_CPRef.Dot(NormCPRefMom);	

	DisCPClu = sqrt(Diff_Clu_CPRef.Mag()*Diff_Clu_CPRef.Mag() - ProDis*ProDis);

	return DisCPClu; 
}

