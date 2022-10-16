/*******************************************************************************
	Stokes.c															6/20/09
	Alexander S Rattner

Defines the basic elements needed to solve 2D Stokes flow
*******************************************************************************/

#include "Stokes.h"

//plugs in the backward diff coeffs into waiting arrays
void GetBDFcoeff(int iOrder, Coeffs *Co)
{
	size_t SZ;
	
	//make space for elements of C
	SZ = 	(size_t)iOrder * sizeof(double);
	Co->dX = (double*) malloc(SZ);
	Co->dRHS = (double*) malloc(SZ);
		
	switch(iOrder)
	{
		case 1:
			Co->dX[0] = 1.0;
			Co->dRHS[0] = 1.0;
			Co->dBDF = 1.0; 
			break;
		case 2:
			Co->dX[0] = -1.0; Co->dX[1] =  2.0;
			Co->dRHS[0] = -0.5; Co->dRHS[1] =  2.0;
			Co->dBDF = 1.5; 
			break;
		case 3:
			Co->dX[0] =  1; Co->dX[1] = -3; Co->dX[2] =  3;
			Co->dRHS[0] =  1/3; Co->dRHS[1] = -1.5; Co->dRHS[2] =  3.0;
			Co->dBDF = 11/6; 
			break;
		case 4:
			Co->dX[0] = -1; Co->dX[1] =  4; Co->dX[2] = -6; Co->dX[3] =  4;
			Co->dRHS[0] = -1/4; Co->dRHS[1] =  4/3; Co->dRHS[2] = -3.0; Co->dRHS[3] =  4.0;
			Co->dBDF = 25/12;
			break;
		default:
			printf("Bad order! Defaulting to order = 1. \n");
			Co->dX[0] = 1.0;
			Co->dRHS[0] = 1.0;
			Co->dBDF = 1.0;
			break;
	}
}

void PrecondL(int N, doublereal* b, doublereal* x, double* dLambda)
{
	fftw_plan pF, pB;
	fftw_complex *cX, *cF;
	size_t SZ;
	int i;
	
	//make room for cX and cF
	SZ = (size_t)N * sizeof(fftw_complex);
	cX = (fftw_complex*) fftw_malloc(SZ);
	cF = (fftw_complex*) fftw_malloc(SZ);	
	
	//make plans
	pF = fftw_plan_dft_1d(N, cX, cF, FFTW_FORWARD, FFTW_MEASURE | FFTW_PRESERVE_INPUT);
	pB = fftw_plan_dft_1d(N, cF, cX, FFTW_BACKWARD, FFTW_MEASURE | FFTW_PRESERVE_INPUT);
	
	//Copy in x
	for (i=0; i<N; i++)
	{
		cX[i][0] = (double) b[i];
		cX[i][1] = 0;
	}
	
	//run preconditioner
	fftw_execute(pF);
	for (i=0; i<N; i++)
	{
		cF[i][0] = cF[i][0] / (dLambda[i] * (double)N);
		cF[i][1] = cF[i][1] / (dLambda[i] * (double)N);	
	}
	fftw_execute(pB);	
	
	//copy result into output
	for (i=0; i<N; i++)
		x[i] = (doublereal) cX[i][0];
	
	//release memory
	fftw_free(cX);
	fftw_free(cF);
	fftw_destroy_plan(pF);
	fftw_destroy_plan(pB);
	
}

//This is the version of the preconditioner that is called by GMRES (indirectly)
void PrecondM(int N, doublereal* b, doublereal* x, double* dL)
{
	size_t SZ;
	int i;
	fftw_complex *cX, *cY;
	fftw_complex *cFX, *cFY;
	fftw_plan pFX, pBX;
	fftw_plan pFY, pBY;
	
	if (N >32)
	{
		//make space for complex variables
		SZ = (size_t)N * sizeof(fftw_complex);
		cX = (fftw_complex*) fftw_malloc(SZ);
		cFX = (fftw_complex*) fftw_malloc(SZ);
		cY = (fftw_complex*) fftw_malloc(SZ);
		cFY = (fftw_complex*) fftw_malloc(SZ);
		
		//make plans
		pFX = fftw_plan_dft_1d(N, cX, cFX, FFTW_FORWARD, FFTW_MEASURE | FFTW_PRESERVE_INPUT);
		pBX = fftw_plan_dft_1d(N, cFX, cX, FFTW_BACKWARD, FFTW_MEASURE | FFTW_PRESERVE_INPUT);
		pFY = fftw_plan_dft_1d(N, cY, cFY, FFTW_FORWARD, FFTW_MEASURE | FFTW_PRESERVE_INPUT);
		pBY = fftw_plan_dft_1d(N, cFY, cY, FFTW_BACKWARD, FFTW_MEASURE | FFTW_PRESERVE_INPUT);
		
		//copy in input vector
		for (i=0; i<N; i++)
		{
			cX[i][0] = b[i]; cX[i][1] = 0;
			cY[i][0] = b[N+i]; cY[i][1] = 0;
		}
		
		//run forward transforms
		fftw_execute(pFX);
		fftw_execute(pFY);
		
		//scale by L/N
		for (i=0; i<N; i++)
		{
			cFX[i][0] = dL[i]*cFX[i][0] / (double)N; cFX[i][1] = dL[i]*cFX[i][1] / (double)N;
			cFY[i][0] = dL[i]*cFY[i][0] / (double)N; cFY[i][1] = dL[i]*cFY[i][1] / (double)N;
		}
		
		//run backward transforms
		fftw_execute(pBX);
		fftw_execute(pBY);
		
		//copy result into x
		for (i=0; i<N; i++)
		{
			x[i]   = cX[i][0];
			x[N+i] = cY[i][0];
		}
		
		//clean up
		fftw_free(cX); fftw_free(cY);
		fftw_free(cFX); fftw_free(cFY);
		fftw_destroy_plan(pFX); fftw_destroy_plan(pBX);
		fftw_destroy_plan(pFY); fftw_destroy_plan(pBY);
	}
	else //simplified version
	{
		SZ = (size_t)(2*N) * sizeof(doublereal);
		memcpy(x, b, SZ);
	}
}

//updates the STORE, dRHS, and dSigM
void UpdateSTORERHS(VesCase *VC, Coeffs *Co)
{
	int i;
	int N = VC->N;
	size_t SZ;

	//copy results into store
	SZ = (size_t)N * sizeof(double);
	memcpy(VC->STORE.dX, VC->V.C.dXPts, SZ);
	memcpy(VC->STORE.dY, VC->V.C.dYPts, SZ);
	
	//fill in STORE and dSigM
	for (i=0; i<N; i++)
	{
		VC->dRHSVec[i]   = VC->STORE.dX[i]*Co->dRHS[0];
		VC->dRHSVec[i+N] = VC->STORE.dY[i]*Co->dRHS[0];
		VC->dSigM[i]     = VC->STORE.dSig[i]*Co->dX[0];
	}
}

//This is currently limited to first order
void UpdateFirstM(Ves2D *V, Params *Par, Store *STORE)
{
	int i, jj;
	size_t SZ;
	int iMR, iM, iN, iNV; //I don't really know what these do
	//Time horizon for the first m steps
	double dT = (double)(Par->iOrder-1)*Par->dTs;
	double dTs; // the time step size
	//Difference coefficients
	Coeffs Co;
	//make storage for vesicle at previous time steps
	SimContour *scHist;
	//I'm not sure what these are for
	double *dSigma;
	double dAvgS[] = {0,0,0,0};
	
	//This is some obscure bookkeeping stuff
	iMR = (int)ceil( ((double)Par->iM)/32 ) * (int)pow(Par->iOrder, 2);
	iM = iMR*(Par->iOrder-1);
	//This will be changed eventually
	iN = V->C.N; 
	iNV = 1; 
	
	//Some bonus bookkeeping stuff
	int iTOrder; //the lowest order of something
	
	//find initial time step size
	if (Par->iOrder > 1)
	{	dTs = dT / (double)iM; }
	else
	{	dTs = dT;}
	
	//prepare some complex thingy
	SZ = (size_t)iN * sizeof(fftw_complex);
	
	//Make space for vesicle history
	SZ = (size_t)(iM+1) *sizeof(SimContour);
	scHist = (SimContour*) malloc(SZ);
	for (i=0; i<(iM+1); i++)
		scHist[i] = Make_SimContour(iN);
	
	//Copy first vesicle state into history
	SZ = (size_t)iN * sizeof(double);
	memcpy(scHist[0].dXPts, V->C.dXPts, SZ);
	memcpy(scHist[0].dYPts, V->C.dYPts, SZ);
	
	//Initialize more stuff
	dSigma = (double*) calloc((size_t)iN, sizeof(double));
	
	//THE MAIN LOOP!!, this doesn't get called for lower orders
	for (jj = 2; jj<=(iM+1); jj++)
	{
		//find the lowest order of something
		if ((jj-1)<=Par->iOrder)
		{	iTOrder = (jj-1);}
		else
		{	iTOrder = Par->iOrder;}
		
		//get the correct coefficients for iTOrder
		GetBDFcoeff(iTOrder, &Co);
		
		//Some matrix multiplications for differences:
		
		//Clean up the old coefficients
		Dest_Coeffs(&Co);
	}
	
	//compile returns into the Vals struct
	SZ = (size_t)iN*sizeof(double);
	memcpy(STORE->dX, scHist[0].dXPts, SZ);
	memcpy(STORE->dY, scHist[0].dYPts, SZ);
	memcpy(STORE->dSig, dSigma, SZ);
	SZ = (size_t)4*sizeof(double);
	memcpy(STORE->dAvgS, dAvgS, SZ);
	
	//free up used space
	for (i=0; i<(iM+1); i++)
		Dest_SimContour(scHist[i]);
	free(scHist);
	free(dSigma);
}

// Computes the Single Layer integral of Stokes kernel at target points dTX dTY
//from the source points dSX/dSY
void ExactStokes(Force *F, double *dTX, double *dTY, int NT, double *dSX, double *dSY, int NS, Force *FV)
{
	int j, i, k; //k becomes the length of iInd
	double *dDis, *dDiffX, *dDiffY;
	size_t SZ;
	int *iInd;
	
	//make space for dDis, dDiffX, dDiffY
	SZ = (size_t)NS * sizeof(double);
	dDis = (double*) malloc(SZ);
	dDiffX = (double*) malloc(SZ);
	dDiffY = (double*) malloc(SZ);
	//make room for iInd
	SZ = (size_t)NS * sizeof(int);
	iInd = (int*)malloc(SZ);
	
	//Fill in F
	for (j=0; j<NT; j++)
	{
		//fill in dDis and dDiffX/Y
		for (i=0; i<NS; i++)
		{
			dDiffX[i] = dSX[i] - dTX[j];
			dDiffY[i] = dSY[i] - dTY[j];
			dDis[i] = sqrt( dDiffX[i]*dDiffX[i] + dDiffY[i]*dDiffY[i] );
		}
		
		//find nonzeros in dDis
		k = 0;
		for (i=0; i<NS; i++) //k becomes the length of iInd
			if( dDis[i] != 0)
			{	iInd[k] = i; k++; }
		
		//fill in F
		F->dFX[j] = 0; F->dFY[j] = 0; //initialize these to 0
		for (i=0; i<k; i++)
		{
			F->dFX[j] -= FV->dFX[iInd[i]]*log(dDis[iInd[i]]);
			F->dFY[j] -= FV->dFY[iInd[i]]*log(dDis[iInd[i]]);
		
		}
		for (i=0; i<k; i++)
		{
			F->dFX[j] += ( pow(dDiffX[iInd[i]] / dDis[iInd[i]], 2) * FV->dFX[iInd[i]] ) +
			             ( dDiffX[iInd[i]] * dDiffY[iInd[i]] * FV->dFY[iInd[i]] / pow(dDis[iInd[i]], 2)  );
			F->dFY[j] += ( pow(dDiffY[iInd[i]] / dDis[iInd[i]], 2) * FV->dFY[iInd[i]] ) +
						 ( dDiffX[iInd[i]] * dDiffY[iInd[i]] * FV->dFX[iInd[i]] / pow(dDis[iInd[i]], 2)  );
		}
	}
	
	//scaling
	cblas_dscal(NT, 1/(4*M_PI), F->dFX, 1);
	cblas_dscal(NT, 1/(4*M_PI), F->dFY, 1);

	//clean up
	free(dDiffX); free(dDiffY);
	free(dDis);
	free(iInd);
}

//Calculates the interaction force between multiple vesicles
void InteractionForce(int NV, VesCase *VC)
{
	int N;
	int i, j;
	int iTotPoints;
	size_t SZ;
	double *dSa, *dISa; //some kind of derivative norm
	double *dDX1, *dDY1, *dDX4, *dDY4; //some sort of weirdly scaled derivative
	double *dSig, *dSigX, *dSigY;
	double *dXAll, *dYAll; //used for interactions
	Force FVAll; //used for interactions
	Force *FV;
	Ves2D *V;

	//initialize dSigX and dSigY
	SZ = (size_t)(VC[0].N) * sizeof(double);
	dSigX = (double*) malloc(SZ);
	dSigY = (double*) malloc(SZ);
	
	//calculate FV for each vesicle
	for (j=0; j<NV; j++)
	{
		N = VC[j].N;
		SZ = (size_t)N * sizeof(double);
		
		//resize dSigX and dSigY (should be extra cheap for first run)
		dSigX = (double*) realloc(dSigX, SZ);
		dSigY = (double*) realloc(dSigY, SZ);
		
		//link to arc lengths, derivatives, dSig, FV, and V
		dSa  = VC[j].V.C.dSa;
		dISa = VC[j].V.C.dISa;
		dDX1 = VC[j].V.C.dDXISa;
		dDY1 = VC[j].V.C.dDYISa;
		dDX4 = VC[j].V.C.dD4XISa;
		dDY4 = VC[j].V.C.dD4YISa;
		dSig = VC[j].dSigM;
		FV = &(VC[j].FV);
		V  = &(VC[j].V);
		
		//fill in dSig*X
		for (i=0; i<N; i++)
		{
			dSigX[i] = dSig[i]*dDX1[i];
			dSigY[i] = dSig[i]*dDY1[i];
		}
		//take its derivative
		DmFT_Double (N, 1, dSigX, dSigX, 1, dISa);
		DmFT_Double (N, 1, dSigY, dSigY, 1, dISa);
			
		//calculate FV
		for (i=0; i<N; i++)
		{
			FV->dFX[i] = (2*M_PI/(double)N) * dSa[i] * ( -(V->dKappa)*dDX4[i] + dSigX[i] );
			FV->dFY[i] = (2*M_PI/(double)N) * dSa[i] * ( -(V->dKappa)*dDY4[i] + dSigY[i] );
		}
		
		//calculate self force
		ExactStokes( &(VC[j].FS), VC[j].V.C.dXPts, VC[j].V.C.dYPts, VC[j].N, VC[j].V.C.dXPts, VC[j].V.C.dYPts, VC[j].N, FV);
	}
	
	//put all points and FVs together in big fat vectors
	iTotPoints = 0;
	for (j=0; j<NV; j++)
		iTotPoints += VC[j].N;
	SZ = (size_t)iTotPoints * sizeof(double);
	dXAll = (double*) malloc(SZ);
	dYAll = (double*) malloc(SZ);
	Make_Force(iTotPoints, &FVAll);
	i = 0; //position counter
	for (j=0; j<NV; j++)
	{
		N = VC[j].N;
		V = &(VC[j].V);
		FV = &( VC[j].FV );
		SZ = (size_t)N * sizeof(double);
		memcpy(&(dXAll[i]), V->C.dXPts, SZ); 
		memcpy(&(dYAll[i]), V->C.dYPts, SZ);
		memcpy( &(FVAll.dFX[i]), FV->dFX, SZ);
		memcpy( &(FVAll.dFY[i]), FV->dFY, SZ);
		i += N;
	}
	
	//calculate interaction forces
	for (j=0; j<NV; j++)
		ExactStokes( &(VC[j].FO), VC[j].V.C.dXPts, VC[j].V.C.dYPts, VC[j].N, dXAll, dYAll, iTotPoints, &FVAll);
	
	//we don't do ViscCont (something related to it would go here)
	
	//return memory
	free(dSigX); free(dSigY);
	free(dXAll); free(dYAll);
	Dest_Force(&FVAll);
}

//this solves for mu, it ought to be called after GetInvLHS
void EvalDensity(Boundary *B, VesCases *VCs)
{
	int NV = VCs->NV;
	int NVP; //the number of vesical points
	int NC = B->NC;
	int iNLets = B->iNLets;
	int NP = B->iTotalSize;
	int NT = 2*NP + 3*iNLets*(NC-1); //total length
	int N;
	double *dRHS;
	size_t SZ;
	int i, j, ind;
	double *dVXAll, *dVYAll; //all of the vesicle points
	double *dBXAll, *dBYAll;
	Force FVB, FVAll, *FV; //force from vesicles on boundary
	
	//get effect of vesicles on boundary
	//first do some initialization
	NVP = 0;
	for (i=0; i<NV; i++)
		NVP += VCs->VC[i].N;
	SZ = (size_t)NVP * sizeof(double);
	dVXAll = (double*) malloc(SZ);
	dVYAll = (double*) malloc(SZ);
	Make_Force(NVP, &FVAll);
	
	//then copy everything into big arrays/structs
	ind = 0;
	for (i=0; i<NV; i++)
	{
		//some bookkeeping
		N = VCs->VC[i].N;
		SZ = (size_t)N * sizeof(double);
		FV = &( VCs->VC[i].FV );
		
		memcpy( &(dVXAll[ind]), VCs->VC[i].V.C.dXPts, SZ);
		memcpy( &(dVYAll[ind]), VCs->VC[i].V.C.dYPts, SZ);
		memcpy( &(FVAll.dFX[ind]), FV->dFX, SZ);
		memcpy( &(FVAll.dFY[ind]), FV->dFY, SZ);
		
		ind += N;
	}
	
	//make room for dBXAll/dBYALL, and fill them with dXPts and dYPts from the boundary
	SZ = (size_t)NP * sizeof(double);
	dBXAll = (double*) malloc(SZ);
	dBYAll = (double*) malloc(SZ);
	ind = 0;
	for (i=0; i<NC; i++)
	{
		N = B->C[i].N;
		SZ = (size_t)N * sizeof(double);
		memcpy( &(dBXAll[ind]), B->C[i].dXPts, SZ);
		memcpy( &(dBYAll[ind]), B->C[i].dYPts, SZ);
		ind += N;
	}
	
	Make_Force(NP, &FVB);
	
	//calculate effect of vesicles on boundary
	ExactStokes(&FVB, dBYAll, dBXAll, NP, dVXAll, dVYAll, NVP, &FVAll);
	
	
	//make room for dRHS, and fill it in
	dRHS = calloc(NT, sizeof(double));	
	ind = 0;
	for (i=0; i<NC; i++)
	{
		N = B->C[i].N;
		SZ = (size_t)N * sizeof(double);
		for (j=0; j<N; j++)
		{
			dRHS[2*j     + ind] = B->C[i].dU[j] - FVB.dFX[j + (ind/2)];
			dRHS[2*j + 1 + ind] = B->C[i].dV[j] - FVB.dFY[j + (ind/2)];
		}
		ind += 2*N;
	}
	
	//perform matrix mult
	cblas_dgemv(CblasRowMajor, CblasNoTrans, NT, NT, 1.0, B->dInvLHS, NT, dRHS, 1, 0.0, B->dMu, 1);	
	
	//clean up
	free(dRHS);
	free(dVXAll); free(dVYAll);
	free(dBXAll); free(dBYAll);
	Dest_Force(&FVB); Dest_Force(&FVAll);
}

// KernelS returns the integration kernel for the single layer potential over the contour
//The output is an 2Nx2N matrix
void KernelS(Ves2D *V, SingQuad *SQ, double *dG)
{
	size_t SZ;
	int N = V->C.N;
	int i, j;
	int *iInd;
	double *dXind, *dYind;
	double *dXin, *dYin;
	double *dRho, *dBr;
	double *dWLR, *dWX, *dWY, *dWXY;
	double *dLogPart;
	double *dXPart, *dYPart, *dXYPart;
	
	//make room for weird counting array
	SZ = (size_t)N * sizeof(int);
	iInd = (int*) malloc(SZ);
	
	//make space for dXin and dYin
	SZ = (size_t)SQ->iY * sizeof(double);
	dXind = (double*)malloc(SZ); //reordered XPts and YPts
	dYind = (double*)malloc(SZ);
	dXin = (double*)malloc(SZ);
	dYin = (double*)malloc(SZ);
		
	//make room for Rho
	SZ = (size_t)SQ->iY * sizeof(double);
	dRho = (double*)malloc(SZ);
	dBr = (double*)malloc(SZ);
	
	//make room for these intermediate vectors
	dWLR = (double*)malloc(SZ);
	dWX = (double*)malloc(SZ);
	dWY = (double*)malloc(SZ);
	dWXY = (double*)malloc(SZ);
	SZ = (size_t)N * sizeof(double);
	dLogPart = (double*)malloc(SZ);
	dXPart = (double*)malloc(SZ);
	dYPart = (double*)malloc(SZ);
	dXYPart = (double*)malloc(SZ);
		
	for (j=0; j<N; j++)
	{
		//Make this counter
		for (i=0; i<N; i++)
			iInd[i] = (j + i) % N;
		
		//make reordered dX and dY
		for (i=0; i<N; i++)
		{
			dXind[i] = V->C.dXPts[iInd[i]];
			dYind[i] = V->C.dYPts[iInd[i]];
		}
		
		//do some matrix multiplication for dXin/dYin
		cblas_dgemv(CblasRowMajor, CblasNoTrans, SQ->iY, SQ->iX, 1.0, SQ->dA, SQ->iX, dXind, 1, 0.0, dXin, 1);
		cblas_dgemv(CblasRowMajor, CblasNoTrans, SQ->iY, SQ->iX, 1.0, SQ->dA, SQ->iX, dYind, 1, 0.0, dYin, 1);
		
		//fill in Rho, some kind of distance, and its inverse: dBr
		for (i=0; i<SQ->iY; i++)
		{
			dRho[i] = pow(dXin[i]-V->C.dXPts[j], 2) + pow(dYin[i]-V->C.dYPts[j], 2);
			dBr[i] = 1 / dRho[i];
		}
		
		//fill in dWLR, dW, dWY
		for (i=0; i<SQ->iY; i++)
		{
			dWLR[i] = -0.5*SQ->dWt[i]*log(dRho[i]);
			dWX[i] = SQ->dWt[i] * dBr[i] * pow(V->C.dXPts[j]-dXin[i], 2);
			dWY[i] = SQ->dWt[i] * dBr[i] * pow(V->C.dYPts[j]-dYin[i], 2);
			dWXY[i] = SQ->dWt[i] * dBr[i] *(V->C.dXPts[j]-dXin[i])*(V->C.dYPts[j]-dYin[i]);
		}
		
		//Compute LogPart, dXPart, dYPart, dXYPart
		cblas_dgemv(CblasRowMajor, CblasTrans, SQ->iY, SQ->iX, 1.0, SQ->dA, SQ->iX, dWLR, 1, 0.0, dLogPart, 1);
		cblas_dgemv(CblasRowMajor, CblasTrans, SQ->iY, SQ->iX, 1.0, SQ->dA, SQ->iX, dWX, 1, 0.0, dXPart, 1);
		cblas_dgemv(CblasRowMajor, CblasTrans, SQ->iY, SQ->iX, 1.0, SQ->dA, SQ->iX, dWY, 1, 0.0, dYPart, 1);
		cblas_dgemv(CblasRowMajor, CblasTrans, SQ->iY, SQ->iX, 1.0, SQ->dA, SQ->iX, dWXY, 1, 0.0, dXYPart, 1);
		
		//Fill in G
		for (i=0; i<N; i++)
		{
			dG[2*N*j     +   iInd[i]] = V->C.dSa[iInd[i]] * (dLogPart[i] + dXPart[i]);
			dG[2*N*j     + N+iInd[i]] = V->C.dSa[iInd[i]] * dXYPart[i];
			dG[2*N*(j+N) +   iInd[i]] = V->C.dSa[iInd[i]] * dXYPart[i];
			dG[2*N*(j+N) + N+iInd[i]] = V->C.dSa[iInd[i]] * (dLogPart[i] + dYPart[i]);
		}
	}
	
	//clear memory
	free(iInd);
	free(dXind); free(dYind);
	free(dXin); free(dYin);
	free(dRho); free(dBr);
	free(dWLR); free(dWX); free(dWY); free(dWXY);
	free(dLogPart);
	free(dXPart); free(dYPart); free(dXYPart);
}

//double layer kernel
void KernelD(Ves2D *V, double *dD)
{
	int N = V->C.N;
	int i, j;
	size_t SZ;
	double *dRX, *dRY;
	double *dRho;
	double *dCo;
	
	//make space for some vars
	SZ = (size_t)N * sizeof(double);
	dRX = (double*) malloc(SZ);
	dRY = (double*) malloc(SZ);
	dRho = (double*) malloc(SZ);
	dCo = (double*) malloc(SZ);
	
	//Fiil in D
	for (j=0; j<N; j++)
	{
		//calculate distances and differences
		for (i=0; i<N; i++)
		{
			dRX[i] = V->C.dXPts[j] - V->C.dXPts[i];
			dRY[i] = V->C.dYPts[j] - V->C.dYPts[i];
			if (i != j)
			{	dRho[i] = pow(dRX[i], 2) + pow(dRY[i], 2); }
			else //i==j
			{	dRho[i] = 1; }
		}
		
		//calculate coefficients, are these units correct?
		for (i=0; i<N; i++)
			dCo[i] = ( V->C.dSa[i] * (dRX[i]*V->C.dNormX[i] + dRY[i]*V->C.dNormY[i]) ) /
			          ( pow(dRho[i], 2) * M_PI );
		
		//Fill in D
		for (i=0; i<N; i++)
		{
			if (i != j)
			{
				dD[2*N*j     +   i] = dCo[i]*dRX[i]*dRX[i] * ((2*M_PI)/N);
				dD[2*N*j     + N+i] = dCo[i]*dRX[i]*dRY[i] * ((2*M_PI)/N);
				dD[2*N*(j+N) +   i] = dCo[i]*dRY[i]*dRX[i] * ((2*M_PI)/N);
				dD[2*N*(j+N) + N+i] = dCo[i]*dRY[i]*dRY[i] * ((2*M_PI)/N);
			}
			else
			{
				dD[2*N*j     +   j] = ( -V->C.dSa[j]*V->C.dKap[j] * V->C.dTangX[j]*V->C.dTangX[j] )/ N;
				dD[2*N*j     + N+j] = ( -V->C.dSa[j]*V->C.dKap[j] * V->C.dTangX[j]*V->C.dTangY[j] )/ N;
				dD[2*N*(j+N) +   j] = ( -V->C.dSa[j]*V->C.dKap[j] * V->C.dTangY[j]*V->C.dTangX[j] )/ N;
				dD[2*N*(j+N) + N+j] = ( -V->C.dSa[j]*V->C.dKap[j] * V->C.dTangY[j]*V->C.dTangY[j] )/ N;
			}
		}
	}
	
	//free up space
	free(dRX); free(dRY);
	free(dRho);
	free(dCo);
}
void TimeStepper(Params *Par, Coeffs *Co, Force *FFar, double *dG, double *dD, Ves2D *V, double *dSig, 
                 Force *FOut, double *dRHSVec, GMRESParams *GMParLHS, GMRESParams *GMParTime)
{
	int N = V->C.N;
	int i, j;
	size_t SZ;
	double *dImD;
	double *dQ;
	double *dGVec;
	double *dRHSTension, *dRHSTenSD;
	double *dSigKe;
	double *dRHSX, *dRHSY, *dRHSXY, *dRHSG;
	double *dFs_expX, *dFs_expY;
	double *dXn;
	double *dY1; //used for getting dSig and FOut
	//make space for dImD (eye minus D)
	SZ = (size_t)((2*N)*(2*N)) * sizeof(double);
	dImD = (double*) malloc(SZ);

	//fill in dImD
	for (i=0; i<2*N; i++)
		for (j=0; j<2*N; j++)
		{
			if (i != j) //off-diagonal
			{	dImD[2*N*i + j] = - dD[2*N*i + j]; }
			else //on diagonal
			{	dImD[2*N*i + j] = 1 - dD[2*N*i + j]; }
		}
	
	//make room for dQ
	SZ = (size_t)(2*N) * sizeof(double);
	dQ = (double*) malloc(SZ);
	//do matrix mult for dQ
	cblas_dgemv(CblasRowMajor, CblasNoTrans, 2*N, 2*N, 1.0, dImD, 2*N, dRHSVec, 1, 0.0, dQ, 1);
	//add in Far Field Force
	for (i=0; i<N; i++)
	{
		dQ[i]   = dQ[i]   + FFar->dFX[i];
		dQ[i+N] = dQ[i+N] + FFar->dFY[i];
	}
	
	//Fill in dGVec
	dGVec = SurfaceDiv(N, dRHSVec, V->C.dDXISa, V->C.dDYISa, V->C.dISa);
	
	//No entropic case for now
	
	//For now only the implicit SchurComp problem type
		
	//fill in RHSTension
	SZ = (size_t)(2*N) * sizeof(double);
	dRHSTension = (double*) malloc(SZ);
	memcpy(dRHSTension, dQ, SZ);
	dRHSTenSD = SurfaceDiv(N, dRHSTension, V->C.dDXISa, V->C.dDYISa, V->C.dISa);
	//sorry about the realloc, it's there in the matlab code
	SZ = (size_t)N * sizeof(double);
	dRHSTension = (double*)realloc(dRHSTension, SZ);
	for (i=0; i<N; i++)
		dRHSTension[i] = dGVec[i] - dRHSTenSD[i];
	
	//make room for dSigKe, and solve it
	dSigKe = (double*) calloc(N, sizeof(double));
	
	//Solve for dSigKe
	RunGMRES(N, dSigKe, dRHSTension, GMParLHS);
	
	//allocate RHS
	SZ = (size_t)N * sizeof(double);
	dRHSX = (double*) malloc(SZ);
	dRHSY = (double*) malloc(SZ);
	dRHSXY = (double*) malloc(2*SZ);
	dRHSG = (double*) malloc(2*SZ);
	
	//fill in RHS
	for (i=0; i<N; i++)
	{
		dRHSX[i] = dSigKe[i]*V->C.dDXISa[i];
		dRHSY[i] = dSigKe[i]*V->C.dDYISa[i];
	}
	
	//Sigke gets scaled here
	cblas_dscal(N, 1/Par->dTs, dSigKe, 1);
	
	//run DmFTs
	DmFT_Double (N, 1, dRHSX, dRHSX, 1, V->C.dISa);
	DmFT_Double (N, 1, dRHSY, dRHSY, 1, V->C.dISa);
	
	//fill in dRHSG, and run the matrix multiply and vector add
	memcpy(&(dRHSXY[0]), dRHSX, SZ);
	memcpy(&(dRHSXY[N]), dRHSY, SZ);
	cblas_dgemv(CblasRowMajor, CblasNoTrans, 2*N, 2*N, 1.0, dG, 2*N, dRHSXY, 1, 0.0, dRHSG, 1);
	for (i=0; i<(2*N); i++)
		dRHSG[i] = dRHSG[i] + dQ[i];
	
	//No entropic forcing for now
	
	//I'm not quite sure what these are
	//but they are just a scaled version of dRHSX/Y
	dFs_expX = (double*) malloc(SZ);
	dFs_expY = (double*) malloc(SZ);
	memcpy(dFs_expX, dRHSX, SZ);
	memcpy(dFs_expY, dRHSY, SZ);
	cblas_dscal(N, 1/Par->dTs, dFs_expX, 1);
	cblas_dscal(N, 1/Par->dTs, dFs_expY, 1);

	//run major GMRES solve
	SZ = (size_t)N * sizeof(double);
	dXn = (double*) malloc(2*SZ);
	//copy in current position to initialize
	memcpy(&(dXn[0]), V->C.dXPts, SZ);
	memcpy(&(dXn[N]), V->C.dYPts, SZ);

	//Major solve
	RunGMRES((2*N), dXn, dRHSG, GMParTime);
	
	//fill in dSig and FOut from TimeMatVec
	//dY1
	SZ = (size_t)(2*N) * sizeof(double);
	dY1 = (double*) malloc(SZ);
	TimeMatVec(dXn, V, dG, dD, Par, Co, dY1, dSig, FOut->dFX, FOut->dFY, GMParLHS);
	
	//Copy GMRES result into Ves2D
	SZ = (size_t)N * sizeof(double);
	memcpy(V->C.dXPts, &(dXn[0]), SZ);
	memcpy(V->C.dYPts, &(dXn[N]), SZ);

	//sum up forces and change direction, sum up sigmas
	for (i=0; i<N; i++)
	{
		FOut->dFX[i] = -(FOut->dFX[i]+dFs_expX[i]);
		FOut->dFY[i] = -(FOut->dFY[i]+dFs_expY[i]);
		dSig[i] = dSig[i]+dSigKe[i];
	}
	
	//return memory
	free(dImD);
	free(dGVec);
	free(dQ);
	free(dRHSTension); free(dRHSTenSD);
	free(dSigKe);
	free(dRHSX); free(dRHSY); free(dRHSXY); free(dRHSG);
	free(dFs_expX); free(dFs_expY);
	free(dXn);
	free(dY1);
}

//I'm not really sure what this is all about
//N is the size of dTangX, which is half of the size of dX
double *SurfaceDiv(int N, double *dX, double *dTangX, double *dTangY, double *dS)
{
	int i;
	double *dX1, *dX2;
	size_t SZ;
	double *dOut;
	
	//make some space
	SZ = (size_t)N * sizeof(double);
	dX1 = (double*) malloc(SZ);
	dX2 = (double*) malloc(SZ);
	//copy in two halves of dX
	memcpy(dX1, &(dX[0]), SZ);
	memcpy(dX2, &(dX[N]), SZ);
	
	//Run DmFTs
	DmFT_Double (N, 1, dX1, dX1, 1, dS);
	DmFT_Double (N, 1, dX2, dX2, 1, dS);
	
	//make space for dOut and fill it in
	SZ = (size_t)N * sizeof(double);
	dOut = (double*) malloc(SZ);
	
	for (i=0; i<N; i++)
		dOut[i] = dX1[i]*dTangX[i] + dX2[i]*dTangY[i];
	
	//free up some memory
	free(dX1); free(dX2);
	return(dOut);
}

void LHSMatVec(double* dSigma, Ves2D *V, double *dG, double *dD, double *dOut)
{
	int i;
	int N = V->C.N;
	double *dSigX, *dSigY, *dDSigXY;
	double *dGfs;
	double *dTemp;
	size_t SZ;
	
	//make space for the sigXY vectors and dGfs
	SZ = (size_t)N * sizeof(double);
	dSigX = (double*) malloc(SZ);
	dSigY = (double*) malloc(SZ);
	SZ = (size_t)(2*N) * sizeof(double);
	dDSigXY = (double*) malloc(SZ);
	dGfs = (double*) malloc(SZ);
	
	//fill in dSigX dSigY, this is equivalent to the matrix operation
	for (i=0; i<N; i++)
	{
		dSigX[i] = dSigma[i]*V->C.dDXISa[i];
		dSigY[i] = dSigma[i]*V->C.dDYISa[i];
	}
	
	//Run first derivatives on those two
	DmFT_Double (N, 1, dSigX, dSigX, 1, V->C.dISa);
	DmFT_Double (N, 1, dSigY, dSigY, 1, V->C.dISa);
	//copy them into dDSigXY
	SZ = (size_t)N * sizeof(double);
	memcpy(&(dDSigXY[0]), dSigX, SZ);
	memcpy(&(dDSigXY[N]), dSigY, SZ);
	
	//multiply for dGfs;
	cblas_dgemv(CblasRowMajor, CblasNoTrans, 2*N, 2*N, 1.0, dG, 2*N, dDSigXY, 1, 0.0, dGfs, 1);
	
	//run surface div
	dTemp = SurfaceDiv(N, dGfs, V->C.dDXISa, V->C.dDYISa, V->C.dISa);
	
	//copy memory into output array
	memcpy(dOut, dTemp, SZ);
	
	//free up memory
	free(dSigX); free(dSigY);
	free(dDSigXY);
	free(dTemp);
	free(dGfs);
}

//TimeMatVec computes the single layer integral over X. To do this, it also calculates the tension.
//Right now there is only support for SchurComp
void TimeMatVec(double *dX, Ves2D *V, double* dG, double *dD, Params *Par, Coeffs *Co, double *dVal, double *dSig, double *dFX, double *dFY, GMRESParams *GMParLHS)
{
	int N = V->C.N;
	int i;
	double *dFkX, *dFkY, *dFkXY;
	double *dXX, *dXY;
	double *dVal1, *dVal1_1, *dVal1_2;
	double *dVal1X, *dVal1Y; //two halves of dVal1
	double *dVal2, *dVal2X, *dVal2Y;
	double *dVal2XY, *dVal2XYG;
	size_t SZ;
	
	if (strcmp(Par->sCase, "schurComp") == 0)
	{
		//copy dX into dXX and dXY
		SZ = (size_t)N * sizeof(double);
		dXX = (double*) malloc(SZ);
		dXY = (double*) malloc(SZ);
		memcpy(dXX, &(dX[0]), SZ);
		memcpy(dXY, &(dX[N]), SZ);
		
		//allocate and fill in dFkX/Y, then scale them
		dFkX = (double*) malloc(SZ);
		dFkY = (double*) malloc(SZ);
		dFkXY = (double*) malloc(2*SZ);
		DmFT_Double (N, 4, dXX, dFkX, 1, V->C.dISa);
		DmFT_Double (N, 4, dXY, dFkY, 1, V->C.dISa);
		cblas_dscal(N, -V->dKappa, dFkX, 1);
		cblas_dscal(N, -V->dKappa, dFkY, 1);
		
		//copy into dFkXY
		memcpy(&(dFkXY[0]), dFkX, SZ);
		memcpy(&(dFkXY[N]), dFkY, SZ);
		
		//make room for dVal1s
		SZ = (size_t)(2*N) * sizeof(double);
		dVal1 = (double*) malloc(SZ);
		dVal1_1 = (double*) malloc(SZ);
		dVal1_2 = (double*) malloc(SZ);
		
		//fill in dVal1
		cblas_dgemv(CblasRowMajor, CblasNoTrans, 2*N, 2*N, 1.0, dD, 2*N, dX, 1, 0.0, dVal1_1, 1);
		cblas_dscal(2*N, Co->dBDF, dVal1_1, 1);
		cblas_dgemv(CblasRowMajor, CblasNoTrans, 2*N, 2*N, 1.0, dG, 2*N, dFkXY, 1, 0.0, dVal1_2, 1);
		cblas_dscal(2*N, Par->dTs, dVal1_2, 1);
		for (i=0; i<(2*N); i++)
			dVal1[i] = dVal1_1[i] + dVal1_2[i];
		
		//fill in dVal1X, dVal1Y (copies of dVal1)
		SZ = (size_t)N * sizeof(double);
		dVal1X = (double*) malloc(SZ);
		dVal1Y = (double*) malloc(SZ);
		memcpy(dVal1X, &(dVal1[0]), SZ);
		memcpy(dVal1Y, &(dVal1[N]), SZ);
		
		//take their first derivatives
		DmFT_Double (N, 1, dVal1X, dVal1X, 1, V->C.dISa);
		DmFT_Double (N, 1, dVal1Y, dVal1Y, 1, V->C.dISa);
		
		//fill in dVal2
		dVal2 = malloc(SZ);
		for (i=0; i<N; i++)
			dVal2[i] = V->C.dDXISa[i]*dVal1X[i] + V->C.dDYISa[i]*dVal1Y[i];
		
		//solve for dSig
		RunGMRES(N, dSig, dVal2, GMParLHS);
		
		//dVal2 is used for another purpose now
		
		//fill in dVal2X/Y with dX1
		dVal2X = (double*) malloc(SZ);
		dVal2Y = (double*) malloc(SZ);
		memcpy(dVal2X, V->C.dDXISa, SZ);
		memcpy(dVal2Y, V->C.dDYISa, SZ);
		
		//scale by sig
		for (i=0; i<N; i++)
		{
			dVal2X[i] = dVal2X[i]*dSig[i];
			dVal2Y[i] = dVal2Y[i]*dSig[i];
		}
		
		//run first derivatives
		DmFT_Double (N, 1, dVal2X, dVal2X, 1, V->C.dISa);
		DmFT_Double (N, 1, dVal2Y, dVal2Y, 1, V->C.dISa);
		
		//fill in net forces: dFX/Y
		for (i=0; i<N; i++)
		{
			dFX[i] = dFkX[i] - (dVal2X[i]/Par->dTs);
			dFY[i] = dFkY[i] - (dVal2Y[i]/Par->dTs);
		}
		
		//run matrix mult on dVal2
		SZ = (size_t)N * sizeof(double);
		dVal2XYG = (double*) malloc(2*SZ);
		dVal2XY = (double*) malloc(2*SZ);
		memcpy(&(dVal2XY[0]), dVal2X, SZ);
		memcpy(&(dVal2XY[N]), dVal2Y, SZ);
		cblas_dgemv(CblasRowMajor, CblasNoTrans, 2*N, 2*N, 1.0, dG, 2*N, dVal2XY, 1, 0.0, dVal2XYG, 1);
		
		//add up dVal
		for (i=0; i<(2*N); i++)
			dVal[i] = Co->dBDF*dX[i] +dVal2XYG[i] - dVal1[i];
		
		//rescale dSig
		cblas_dscal(N, (-1/Par->dTs), dSig, 1);
		
		//return memory
		free(dFkX); free(dFkY); free(dFkXY);
		free(dXX); free(dXY);
		free(dVal1); free(dVal1_1); free(dVal1_2);
		free(dVal1X); free(dVal1Y);
		free(dVal2); free(dVal2X); free(dVal2Y);
		free(dVal2XY); free(dVal2XYG);
	}
}

//Wrapper for the MatVec function (I get called by GMRES)
int LHSMatVecGM(doublereal *alpha, doublereal *x, doublereal *beta, doublereal *y, void** vMData)
{
	//persistent guys
	Ves2D *V;
	double *dG;
	double *dD;
	int N;
	double *dY1;
	size_t SZ;
	
	//read in pack data
	V = (Ves2D*) (vMData[0]);
	dG = (double*) (vMData[1]);
	dD = (double*) (vMData[2]);
	N = V->C.N;

	//make space for two halves of y
	SZ = (size_t)N * sizeof(double);
	dY1 = (double*) malloc(SZ);
		
	//read in data pack
	V = (Ves2D*)(vMData[0]);
	N = V->C.N;
	dG = (double*) (vMData[1]);
	dD = (double*) (vMData[2]);
		
	//compute matvec half
	LHSMatVec(x, V, dG, dD, dY1);
		
	//add up to get result
	cblas_dscal(N, *beta, y, 1);
	cblas_daxpy(N, *alpha, dY1, 1, y, 1);
		
	//free up some space
	free(dY1);
	
	return 0;
}

//Wrapper for TimeMatVec function
int TimeMatVecGM(doublereal *alpha, doublereal *x, doublereal *beta, doublereal *y, void** vMData)
{
	//persistent guys
	Ves2D *V;
	double *dG;
	double *dD;
	Params *Par;
	Coeffs *Co;
	GMRESParams *GMParLHS;
	int N;

	//These are extra results that we only need for the setting case
	size_t SZ;
	double *dSig, *dFX, *dFY;
	//This is the result we do want
	double *dY1;
	
	//Read in data pack
	V = (Ves2D*) (vMData[0]);
	dG = (double*) (vMData[1]);
	dD = (double*) (vMData[2]);
	Par = (Params*) (vMData[3]);
	Co = (Coeffs*) (vMData[4]);
	GMParLHS = (GMRESParams*) (vMData[5]);
	N = V->C.N;

	//make room for extra results
	dSig = (double*) calloc(N, sizeof(double));
	SZ = (size_t)N * sizeof(double);
	dFX = (double*) malloc(SZ);
	dFY = (double*) malloc(SZ);
	dY1 = (double*) malloc(2*SZ);

	//run TimeMatVec
	TimeMatVec(x, V, dG, dD, Par, Co, dY1, dSig, dFX, dFY, GMParLHS);

	//fill in output
	cblas_dscal(2*N, *beta, y, 1);
	cblas_daxpy(2*N, *alpha, dY1, 1, y, 1);
		
	//return some memory		
	free(dSig);
	free(dFX); free(dFY);
	free(dY1);
	
	return 0;
}

//Wrapper for the Preconditioner (I get called by GMRES)
int PrecondLGM(doublereal *x, doublereal *b, void** vPData)
{
	int N;
	double *dLambda;
	
	//get data out of pack
	N  = *( (int*)(vPData[0]) );
	dLambda = (double*)(vPData[1]);
	
	//precond
	PrecondL(N, b, x, dLambda);
	
	return 0;
}

//Wrapper for the M preconditioner (called directly by GMRES)
int PrecondMGM(doublereal *x, doublereal *b, void** vPData)
{
	int N;
	double *dL;
	
	//get data out of pack
	N  = *( (int*)(vPData[0]) );
	dL = (double*)(vPData[1]);
	
	//precond
	PrecondM(N, b, x, dL);
	
	return 0;
}

//Simplified caller for GMRES
void RunGMRES(int iNin, double *dX, double *dB, GMRESParams *GMPar)
{
	integer N = (integer)iNin;
	size_t SZ;
	doublereal *work; //some workspace array
	doublereal *h; //hessenberg matrix
	
	//make space for working arrays
	SZ = (size_t)(GMPar->iLdw * (GMPar->iRestart+4)) * sizeof(doublereal);
	work = (doublereal*) malloc(SZ);
	SZ = (size_t)(GMPar->iLdh * (GMPar->iRestart+2)) * sizeof(doublereal);
	h = (doublereal*) malloc(SZ);
	
	//reset these limits (they store the final iteration count and final residual)
	GMPar->iIterAct = GMPar->iIterMax;
	GMPar->dResAct = GMPar->dResMax;
	
	//run solver	
	gmres(&N, (doublereal*)dB, (doublereal*)dX, &(GMPar->iRestart), work, &(GMPar->iLdw), h, 
	      &(GMPar->iLdh), &(GMPar->iIterAct), &(GMPar->dResAct), GMPar->matvec, GMPar->psolve, GMPar->vPData, GMPar->vMData, &(GMPar->iInfo));
	
	//clean up
	free(work);
	free(h);
}

//calculate average stress, right now this only supports single vesicles
void AvgStress(Ves2D *V, Force *FTJ, double *dAvgS)
{
	int N = V->C.N;
	int i;
	
	//clear dAvgS (always a 4 long array)
	dAvgS[0] = 0; dAvgS[1] = 0; dAvgS[2] = 0; dAvgS[3] = 0;
	
	//Add up stresses
	for (i=0; i<N; i++)
	{
		dAvgS[0] = dAvgS[0] + FTJ->dFX[i]*V->C.dXPts[i]*V->C.dSa[i];
		dAvgS[1] = dAvgS[1] + FTJ->dFY[i]*V->C.dXPts[i]*V->C.dSa[i];
		dAvgS[2] = dAvgS[2] + FTJ->dFX[i]*V->C.dYPts[i]*V->C.dSa[i];
		dAvgS[3] = dAvgS[3] + FTJ->dFY[i]*V->C.dYPts[i]*V->C.dSa[i];
	}
	
	//scale stresses
	cblas_dscal(4, (2*M_PI)/((double)N), dAvgS, 1);
	
	//right now we assume area is 1
}

//TOP LEVEL SOLVER
//This is currently limited to 1st order flow without viscosity contrast, and I'm pretty sure it's incompressible
void SOLVE(VesCases *VCs, Params *Par, Boundary *B, EvalVel EV)
{
	int	i, j;
	int NV = VCs->NV;
	double dTime = 0;
	Coeffs Co;
	
	//get me some coefficients
	GetBDFcoeff(Par->iOrder, &Co);
	
	//some initialization for all of the vesicles
	for (i=0; i<NV; i++)
	{
		//Get first step data
		UpdateFirstM( &(VCs->VC[i].V), Par, &(VCs->VC[i].STORE) );
		//calculate some geometry properties of the vesicle
		//This is called towards the end of the main loop, but has to be called before the first run
		CurveProps( &(VCs->VC[i].V.C) );
		//fill in STORE and dSigM
		UpdateSTORERHS( &(VCs->VC[i]), &Co);
	}
	
	//PUT THE MAIN LOOP IN HERE
	for (j=Par->iOrder+1; j <=(Par->iM+1); j++)
	{
		dTime = dTime + Par->dTs; //time step
		
		//Calculate self force and interaction forces
		InteractionForce(NV, VCs->VC);
		
		//Solve for mu, and velocities (if the flow is bounded)
		if (B != NULL)
		{
			EvalDensity(B, VCs);
			UpsampleMu(B);
		}
		
		//Calculate the velocities
		for (i=0; i<NV; i++)
			EV(VCs->VC[i].N, VCs->VC[i].V.C.dXPts, VCs->VC[i].V.C.dYPts, VCs->VC[i].V.C.dU, VCs->VC[i].V.C.dV, 0);
		
		//Solve for individual vesicles
		for (i=0; i<NV; i++)		
			SOLVE_VesStep( &(VCs->VC[i]), Par, EV, &Co);
		
		//Write output from this step
		WriteOutput("Step", 2, VCs, dTime);
	}
	
	//Clear up memory
	Dest_Coeffs(&Co);
}

void SOLVE_VesStep(VesCase *VC, Params *Par, EvalVel EV, Coeffs *Co)
{
	int i;
	int N = VC->V.C.N;
	
	//Preconditioner data packs
	VC->GMParLHS.vPData[0] = &N;
	VC->GMParLHS.vPData[1] = VC->dLambda;
	VC->GMParTime.vPData[0] = &N;
	VC->GMParTime.vPData[1] = VC->dL;

	
	//Get far field velocity force
	for (i=0; i<N; i++)
	{
		VC->FFar.dFX[i] = ( 2*Par->dTs*( VC->FO.dFX[i] - VC->FS.dFX[i] + VC->V.C.dU[i] ) )/( 1 + VC->V.dViscCont );
		VC->FFar.dFY[i] = ( 2*Par->dTs*( VC->FO.dFY[i] - VC->FS.dFY[i] + VC->V.C.dV[i] ) )/( 1 + VC->V.dViscCont );
	}
	
	//get kernels
	KernelS( &(VC->V), &(VC->SQ), VC->dG);
	KernelD( &(VC->V), VC->dD);
	
	//scale the kernels
	cblas_dscal( ((2*N)*(2*N)), (2/(1+VC->V.dViscCont)), VC->dG, 1);
	cblas_dscal( ((2*N)*(2*N)), (2*(1-VC->V.dViscCont))/(1+VC->V.dViscCont), VC->dD, 1);
	
	//Matvec data packs
	VC->GMParLHS.vMData[0] = &(VC->V);
	VC->GMParLHS.vMData[1] = VC->dG;
	VC->GMParLHS.vMData[2] = VC->dD;
	VC->GMParTime.vMData[0] = &(VC->V);
	VC->GMParTime.vMData[1] = VC->dG;
	VC->GMParTime.vMData[2] = VC->dD;
	VC->GMParTime.vMData[3] = Par;
	VC->GMParTime.vMData[4] = Co;
	VC->GMParTime.vMData[5] = &(VC->GMParLHS);

	//Integration Step, FOut is the traction jump for single vesicle systems
	TimeStepper(Par, Co, &(VC->FFar), VC->dG, VC->dD, &(VC->V), VC->STORE.dSig, &(VC->FOut), VC->dRHSVec, &(VC->GMParLHS), &(VC->GMParTime) ); 
	
	//update store, RHS, dSigM
	UpdateSTORERHS(VC, Co);
	
	//calculate some geometry properties of the vesicle
	CurveProps( &(VC->V.C) );
}
