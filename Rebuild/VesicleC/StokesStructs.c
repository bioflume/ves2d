/*******************************************************************************
	StokesStructs.c														6/27/09
	Alexander S Rattner

Defines all of the structs needed in Stokes (Stokes is getting bloated) and 
their make/dest functions
*******************************************************************************/

#include "StokesStructs.h"


//MAKERS GO UP HERE

//This just runs the Mallocs needed for the Force struct
void Make_Force(int N, Force *F)
{
	size_t SZ;
	F->N = N;
	
	SZ = (size_t)N * sizeof(double);
	F->dFX = malloc(SZ);
	F->dFY = malloc(SZ);
}

//Initializes a STORE
void Make_Store(int N, Store *STORE)
{
	size_t SZ;
	
	SZ = (size_t)N*sizeof(double);
	STORE->dX = (double*) malloc(SZ);
	STORE->dY = (double*) malloc(SZ);
	STORE->dSig = (double*) malloc(SZ);
	SZ = (size_t)4*sizeof(double);
	STORE->dAvgS = (double*) malloc(SZ);
}

//Initializes a Vesicle Problem
void Make_VesCase(int N, char *sDescriptor, VesCase *VC)
{
	size_t SZ;
	
	VC->N = N;
	
	Make_Vesicle(sDescriptor, N, &(VC->V) );
	
	//prepare self force
	Make_Force(N, &(VC->FS) );
	//prepare other interaction force
	Make_Force(N, &(VC->FO) );
	//prepare far field force
	Make_Force(N, &(VC->FFar) );
	//make room for integration (timestepper) ouput
	Make_Force(N, &(VC->FOut) );
	//this is an intermediate force used by InteractionForce
	Make_Force(N, &(VC->FV) );
	//Double layer force
	Make_Force(N, &(VC->FD) );
	
	//make room for the store
	Make_Store(N, &(VC->STORE) );
	
	//make room for the Kernels
	SZ = (size_t)((2*N)*(2*N)) * sizeof(double);
	VC->dG = (double*) malloc(SZ);
	VC->dD = (double*) malloc(SZ);
	
	//make room for dSigM
	VC->dSigM = (double*)calloc(N, sizeof(double));
	
	//Make space for RHSVec
	SZ = (size_t)(2*N)*sizeof(double);
	VC->dRHSVec = (double*)malloc(SZ);
	
	//initialize dLambda and dL
	VC->dLambda = GetLambda(N);
	VC->dL = GetL(N);
	
	//Initialize SingQuad
	QuadLogSingData((N/2), 8, 2*M_PI, &(VC->SQ) ); //I don't know if the 2*pi ever changes
	
	//Initialize the precond and matvec array
	VC->GMParLHS.vPData  = malloc( 2*sizeof(void*) );
	VC->GMParTime.vPData = malloc( 2*sizeof(void*) );
	VC->GMParLHS.vMData  = malloc( 3*sizeof(void*) );
	VC->GMParTime.vMData = malloc( 6*sizeof(void*) );
}

//EVERYTHING ELSE IN THE MIDDLE

//This makes lambda which is a funky counting vector
double *GetLambda(int N)
{
	int i;
	double *dLambda;
	size_t SZ;
	
	//make room for dLambda
	SZ = (size_t)N * sizeof(double);
	dLambda = malloc(SZ);
	
	//apply the weird counting function
	for (i=1; i<N; i++)
	{
		if (i <= (N/2))
		{ dLambda[i] = -0.5*i*M_PI; }
		else
		{ dLambda[i] = -0.5*(N-i)*M_PI; }
	}
	dLambda[0] = dLambda[1];
	
	return dLambda;
}

//L is another counting vector
double *GetL(int N)
{
	int i;
	double *dL;
	size_t SZ;
	
	SZ = (size_t)N * sizeof(double);
	dL = (double*) malloc(SZ);
	for (i=0; i<N; i++)
	{	
		dL[i] = -(N/2) + i;
		dL[i] = abs(pow(dL[i], 3));
		if (dL[i] == 0)
			dL[i] = 1;
		dL[i] = (M_PI/4) * dL[i];
		dL[i] = 1 / dL[i];
	}
	
	return dL;
}

//Gets weights for solving around some sort of singularity?
void GetWeights(int iQ, int iKer, Weights *W)
{
	double **dXp;
	int iLength[3], iPar[3]; //??
	int iY = 25;
	int iX = 2;
	int i;
	size_t SZ;
	
	//I'm not sure if I should leave this as a fixed size
	dXp = Make2D_double(iY, iX);
	
	
	switch(iKer) //Currently I only have the dat file for iKer=2
	{
		case 0: //I don't have this dat file
			iLength[0] = 2; iLength[1] = 4; iLength[2] = 8;
			iPar[0] = 2; iPar[1] = 4; iPar[2] = 7;
			ReadDat("DAT/nodesr.dat", dXp, iY, iX);
			break;
		case 1: //I don't have this dat file
			iLength[0] = 4; iLength[1] = 8; iLength[2] = 16;
			iPar[0] = 3; iPar[1] = 5; iPar[2] = 10;
			ReadDat("DAT/nodes_sqrtx.dat", dXp, iY, iX);
			break;
		case 2: //I do have this dat file
			iLength[0] = 3; iLength[1] = 7; iLength[2] = 15;
			iPar[0] = 2; iPar[1] = 5; iPar[2] = 10;
			ReadDat("DAT/nodes_logx.dat", dXp, iY, iX);
			break;
	}
	
	switch(iQ)
	{
		case 4:
			SZ = (size_t)iLength[0] * sizeof(double);
			W->dV = malloc(SZ);
			W->dU = malloc(SZ);
			for (i=0; i<iLength[0]; i++)
			{ 
				W->dV[i] = dXp[i][0];
				W->dU[i] = dXp[i][1];
			}
			W->iA = iPar[0];
			break;
		case 8:
			SZ = (size_t)iLength[1] * sizeof(double);
			W->dV = malloc(SZ);
			W->dU = malloc(SZ);
			for (i=0; i<(iLength[1]); i++)
			{ 
				W->dV[i] = dXp[iLength[0]+i][0];
				W->dU[i] = dXp[iLength[0]+i][1];
			}
			W->iA = iPar[1];
			break;
		case 16:
			SZ = (size_t)iLength[2] * sizeof(double);
			W->dV = malloc(SZ);
			W->dU = malloc(SZ);
			for (i=0; i<(iLength[2]); i++)
			{ 
				W->dV[i] = dXp[iLength[0]+iLength[1]+i][0];
				W->dU[i] = dXp[iLength[0]+iLength[1]+i][1];
			}
			W->iA = iPar[2];
			break;
	}
	
	Dest2D_double(dXp, iY);
}

//Gets some more weights
void GetRegWeights(int iA, double *dX, double *dW)
{
	int i, iStart, iEnd;
	int iX=2; int iY=35;
	double **dXp;
	
	//Figure out where to start and end
	switch(iA)
	{
		case  2: iStart =  0; iEnd =  1; break;
		case  3: iStart =  2; iEnd =  4; break;
		case  4: iStart =  5; iEnd =  8; break;
		case  5: iStart =  9; iEnd = 14; break;
		case  7: iStart = 15; iEnd = 22; break;
		case 10: iStart = 23; iEnd = 34; break;
	}
	
	//Read in the weight file
	dXp = Make2D_double(iY, iX);
	ReadDat("DAT/nodes_regular.dat", dXp, iY, iX);
	
	//Fill in dX and dW
	for (i=0; i<=(iEnd-iStart); i++)
	{	
		dX[i] = dXp[i+iStart][0];
		dW[i] = dXp[i+iStart][1];
	}
	
	//clear up dXp
	Dest2D_double(dXp, iY);
}

//gives nodes and weights to integrate logarithmic singularity.
void QuadLogSingData(int iM, int iQ, double dL, SingQuad *SQ)
{
	size_t SZ;
	int iKer = 2;
	double dL1 = 0;
	double dL2 = dL/2;
	//The following are for GetWeights
	Weights W;
	//The following are for GetRegWeights
	double *dX, *dW;
	int iLP1, iLP3, iEnd, iEndOld, iLYt;
	double dH = 1/((double)iM);
	//I'm not certain what the following are for
	double *dEvalP1, *dEvalP3;
	double *dYs, *dYt, *dWt;
	int k;
	double *dA, *dB;
	//Used for traversing B
	int *iPos;
	
	//iLP1 is going to be the same size and W.dV and W.dU
	iLP1 = iQ-1;
	SZ = (size_t)iLP1*sizeof(double);
	dEvalP1 = (double*) malloc(SZ);
	
	//Get a bunch of coefficients and points and stuff
	GetWeights(iQ, iKer, &W);
	
	//set dEvalP1
	for (k=0; k<iLP1; k++)
		dEvalP1[k] = W.dV[k]*dH;
	
	//The size of dX and dW depends on iA
	switch(W.iA)
	{
		case  2: iLP3 =  2; break;
		case  3: iLP3 =  3; break;
		case  4: iLP3 =  4; break;
		case  5: iLP3 =  6; break;
		case  7: iLP3 =  8; break;
		case 10: iLP3 = 12; break;
	}
	SZ = (size_t)iLP3*sizeof(double);
	dX = (double*) malloc(SZ);
	dW = (double*) malloc(SZ);
	dEvalP3 = (double*) malloc(SZ);
	GetRegWeights(W.iA, dX, dW);
	
	//set dEvalP3
	for (k=0; k<iLP3; k++)
		dEvalP3[k] = 1 - dX[k]*dH;
	
	//now we can make dYs and dWt, this is MUCH cleaner in MATLAB!!
	iEnd = (2*(iLP1+iLP3))-1;
	SZ = (size_t)(iEnd+1)*sizeof(double);
	dYs = (double*) malloc(SZ);
	dWt = (double*) malloc(SZ);
	for (k=0; k<iLP1; k++)
	{
		dYs[k] = dL1 + (dL2-dL1)*dEvalP1[k];
		dYs[iEnd-k] = dL - dYs[k];
		dWt[k] = (dL2-dL1)*dH*W.dU[k];
		dWt[iEnd-k] = dWt[k];
	}
	for (k=iLP1; k<(iLP1+iLP3); k++)
	{
		dYs[k] = dL1 + (dL2-dL1)*dEvalP3[k-iLP1];
		dYs[iEnd-k] = dL - dYs[k];
		dWt[k] = (dL2-dL1)*dH*dW[k-iLP1];
		dWt[iEnd-k] = dWt[k];
	}
	
	//Make interpolant matrix
	dA = GetInterpMat(dYs, 2*(iLP1+iLP3), 2*iM, dL);
	dH = (dL/2)/((double)iM);
	
	//The MATLAB code that I am working from is ugly, I'm sorry for the realloc
	iLYt = iM-(2*W.iA)+1;
	if (iLYt <= 0) //prevent negative sizes
		iLYt = 0;
	dYt = (double*) malloc((size_t)iLYt*sizeof(double));
	for (k=0; k<iLYt; k++)
		dYt[k] = (double)(W.iA+k)*dH;
	iEndOld = iEnd;
	iEnd = 2*(iLP1+iLP3+iLYt) - 1;
	SZ = (size_t)(iEnd+1)*sizeof(double);
	dWt = (double*) realloc(dWt, SZ);
	for (k=0; k<=iEndOld; k++)
		dWt[k] = dWt[k] / (4*M_PI);
	for (k=(iEndOld+1); k<= iEnd; k++)
		dWt[k] = dH / (4*M_PI);
	
	//now make dB
	iLYt = 2*iLYt;
	if (iLYt <= 0) //prevent negative sizes
		iLYt = 0;
	dB = (double*) calloc(iLYt*2*iM, sizeof(double));
	SZ = (size_t)(iLYt)*sizeof(int);
	iPos = (int*) malloc(SZ);
	
	//Some counting thing
	for (k=0; k<(iLYt/2); k++)
		iPos[k] = W.iA+k;
	for (k=(iLYt/2); k<iLYt; k++)
		iPos[k] = (iM+W.iA) + (k-(iLYt/2));
	
	for (k=0; k<iLYt; k++)
		dB[2*iM*k + iPos[k]] = 1;
	
	//Another realloc to append B to A
	SZ = (size_t)( ( 2*iM )*( 2*(iLP1+iLP3) + iLYt ) ) * sizeof(double); //size of A+B
	dA = (double*) realloc(dA, SZ);
	SZ = (size_t)(2*iM*iLYt)* sizeof(double);//size of B
	memcpy(&(dA[2*(iLP1+iLP3)*2*iM]), dB, SZ);
	
	//clear out mallocs
	Dest_Weights(&W);
	free(dX); free(dW);
	free(dEvalP1); free(dEvalP3);
	free(dYs); free(dYt);
	free(dB); free(iPos);
	
	//piece together return struct
	SQ->dA = dA;
	SQ->dWt = dWt;
	SQ->iY = 2*(iLP1+iLP3) + iLYt;
	SQ->iX = 2*iM;
}

//given a set of x,y points and a velocity field ID, this fills in the velocities
void EvalFFVelocity(int N, double *dX, double *dY, double *dU, double *dV, int iMode, ...)
{
	int i;
	//handle multiple arguments
	va_list vaArgs;
	static VelField *VF;
	
	if (iMode == 1) //initialization case
	{
		va_start (vaArgs, iMode);    
		VF = va_arg(vaArgs, VelField*);
		va_end(vaArgs);
	}
	else //regular velocity calculation
	{
		//apply the particular velocity field
		if ( !strcmp("Shear", VF->sFieldFun) )
		{		
			for (i=0; i<N; i++)
			{
				dU[i] = dY[i]*VF->dMisc[0];
				dV[i] = 0;
			}
		}
		
	}
}


//DESTROYERS GO DOWN HERE
void Dest_SingQuad(SingQuad *SQ)
{
	free(SQ->dA);
	free(SQ->dWt);
}

void Dest_Coeffs(Coeffs *Co)
{
	free(Co->dX);
	free(Co->dRHS);
}

void Dest_Weights(Weights *W)
{
	free(W->dV); free(W->dU);
}

void Dest_Store(Store *STORE)
{
	free(STORE->dX);
	free(STORE->dY);
	free(STORE->dSig);
	free(STORE->dAvgS);
}

void Dest_Force(Force *F)
{
	free(F->dFX);
	free(F->dFY);
}

void Dest_VesCase(VesCase *VC)
{
	Dest_Vesicle( &(VC->V) );
	Dest_Force( &(VC->FS) );
	Dest_Force( &(VC->FO) );
	Dest_Force( &(VC->FFar) );
	Dest_Force( &(VC->FOut) );
	Dest_Force( &(VC->FV) );
	Dest_Force( &(VC->FD) );
	Dest_SingQuad( &(VC->SQ) );
	Dest_Store( &(VC->STORE) );
	free(VC->dG);
	free(VC->dD);
	free(VC->dSigM);
	free(VC->dRHSVec);
	free(VC->dLambda);
	free(VC->dL);
	free( VC->GMParLHS.vPData );
	free( VC->GMParTime.vPData );
	free( VC->GMParLHS.vMData );
	free( VC->GMParTime.vMData );
}

void Dest_VesCases(VesCases *VCs)
{
	int i;
	int NV = VCs->NV;
	for (i=0; i<NV; i++)
		Dest_VesCase( &(VCs->VC[i]) );
	free(VCs->VC);
}

void Dest_VelField(VelField *VF)
{
	free(VF->dMisc);
}
