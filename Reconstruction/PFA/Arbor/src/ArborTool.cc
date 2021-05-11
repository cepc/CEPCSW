#include "ArborTool.h"
#include <TMath.h>

using namespace std;

std::vector<int>SortMeasure( std::vector<float> Measure, int ControlOrderFlag )
{

	std::vector<int> objindex;
	int Nobj = Measure.size();

	for(int k = 0; k < Nobj; k++)
	{
		objindex.push_back(k);
	}

	int FlagSwapOrder = 1;
	float SwapMeasure = 0;
	int SwapIndex = 0;

	for(int i = 0; i < Nobj && FlagSwapOrder; i++)
	{
		FlagSwapOrder = 0;
		for(int j = 0; j < Nobj - 1; j++)
		{
			if((Measure[j] < Measure[j+1] && ControlOrderFlag) || (Measure[j] > Measure[j+1] && !ControlOrderFlag) )
			{
				FlagSwapOrder = 1;
				SwapMeasure = Measure[j];
				Measure[j] = Measure[j+1];
				Measure[j+1] = SwapMeasure;

				SwapIndex = objindex[j];
				objindex[j] = objindex[j+1];
				objindex[j+1] = SwapIndex;
			}
		}
	}

	return objindex;
}

TMatrixF MatrixSummarize( TMatrixF inputMatrix )
{

	int Nrow = inputMatrix.GetNrows();
	int Ncol = inputMatrix.GetNcols();

	TMatrixF tmpMatrix(Nrow, Ncol); 

	for(int i0 = 0; i0 < Nrow; i0 ++)
	{
		for(int j0 = i0; j0 < Ncol; j0 ++)
		{
			//		if( fabs(inputMatrix(i0, j0) - 2) < 0.2 || fabs(inputMatrix(i0, j0) - 10 ) < 0.2 )	//Case 2, 3: Begin-End connector
			if(inputMatrix(i0, j0))
			{
				tmpMatrix(i0, j0) = 1;
				tmpMatrix(j0, i0) = 1;
			}
			else 
			{
				tmpMatrix(i0, j0) = 0;
				tmpMatrix(j0, i0) = 0;
			}
		}
	}

	int PreviousLinks = -1;
	int CurrentLinks = 0;
	int symloopsize = 0;
	vector <int> Indirectlinks;
	int tmpI(0);
	int tmpJ(0);

	// cout<<"Matrix Type: "<<Nrow<<" * "<<Ncol<<endl;

	if( Nrow == Ncol )
	{

		while( CurrentLinks > PreviousLinks )
		{
			PreviousLinks = 0;
			CurrentLinks = 0;

			for(int i = 0; i < Nrow; i ++)
			{
				for(int j = 0; j < Ncol; j ++)
				{
					if( tmpMatrix(i, j) > 0.1 )
						PreviousLinks ++;
				}
			}

			for(int k = 0; k < Nrow; k ++)
			{
				for(int l = 0; l < Ncol; l ++)
				{
					if( tmpMatrix(k, l) > 0.1) Indirectlinks.push_back(l);
				}
				symloopsize = Indirectlinks.size();

				for(int l1(0); l1 < symloopsize; l1 ++)
				{
					tmpI = Indirectlinks[l1];
					for(int m1=l1 + 1; m1 < symloopsize; m1++)
					{
						tmpJ = Indirectlinks[m1];
						tmpMatrix(tmpI, tmpJ) = 1;
						tmpMatrix(tmpJ, tmpI) = 1;
					}
				}
				Indirectlinks.clear();
			}

			for(int u = 0; u < Nrow; u++)
			{
				for(int v = 0; v < Ncol; v++)
				{
					if( tmpMatrix(u, v) > 0.1)
						CurrentLinks ++;
				}
			}

			// cout<<"Link Matrix Symmetrize Loop, PreviousFlag = "<<PreviousLinks<<", CurrentFlag = "<<CurrentLinks<<" of D = "<<Nrow<<" Matrix" <<endl;

		}
	}

	return tmpMatrix; 
}

branchcoll ArborBranchMerge(branchcoll inputbranches, TMatrixF ConnectorMatrix)	//ABM
{
	branchcoll outputbranches; 

	int NinputBranch = inputbranches.size();
	int Nrow = ConnectorMatrix.GetNrows();
	int Ncol = ConnectorMatrix.GetNcols();
	int FlagBranchTouch[Nrow];
	int tmpCellID = 0; 

	if(Ncol != NinputBranch || Nrow != Ncol || Nrow != NinputBranch)
	{
		cout<<"Size of Connector Matrix and inputClusterColl is not match"<<endl;
	}

	for(int i0 = 0; i0 < Nrow; i0++)
	{
		FlagBranchTouch[i0] = 0;
	}

	for(int i1 = 0; i1 < Nrow; i1++)
	{
		if(FlagBranchTouch[i1] == 0)
		{
			branch Mergebranch_A = inputbranches[i1];
			branch a_MergedBranch = Mergebranch_A;
			FlagBranchTouch[i1] = 1;		

			for(int j1 = i1 + 1; j1 < Nrow; j1++)
			{
				if(FlagBranchTouch[j1] == 0 && ConnectorMatrix(i1, j1) > 0.1 )
				{
					branch Mergebranch_B = inputbranches[j1];
					FlagBranchTouch[j1] = 1;
					for(unsigned int k1 = 0; k1 < Mergebranch_B.size(); k1 ++)
					{
						tmpCellID = Mergebranch_B[k1];
						if(find(a_MergedBranch.begin(), a_MergedBranch.end(), tmpCellID) == a_MergedBranch.end())
						{
							a_MergedBranch.push_back(tmpCellID);
						}
					}					
				}
			}
			outputbranches.push_back(a_MergedBranch);
		}
	}	

	return outputbranches; 
}

