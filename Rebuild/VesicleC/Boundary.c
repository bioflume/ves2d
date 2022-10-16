/*******************************************************************************
	Boundary.c															7/24/09
	Alexander S Rattner

Handles all of the fixed boundary specific code.
*******************************************************************************/

#include "Boundary.h"

void Make_Boundary(int iNLets, int NC, Contour *C, Boundary *B)
{
	int i, N, iTotalSize;
	size_t SZ;
	
	//assign number of contours
	B->NC = NC;
	//link in contours
	B->C = C;
	//assign nlets
	B->iNLets = iNLets;
	
	//get curve props for all the contours
	for (i=0; i<NC; i++)
		CurveProps( &(B->C[i]) );
	
	//count up number of points
	iTotalSize = 0;
	for (i=0; i<NC; i++)
		iTotalSize += B->C[i].N;
	B->iTotalSize = iTotalSize;

	//count points in each contour
	SZ = (size_t)NC * sizeof(int);
	B->NP = (int*)malloc(SZ);
	for (i=0; i<NC; i++)
		B->NP[i] = B->C[i].N;
	
	//make room for dInvLHS
	SZ = (size_t)pow((2*iTotalSize+3*iNLets*(NC-1)), 2) * sizeof(double);
	B->dInvLHS = (double*) malloc(SZ);
	
	//make room for dMu
	SZ = (size_t)(2*iTotalSize+3*iNLets*(NC-1)) * sizeof(double);
	B->dMu = (double*) malloc(SZ);
	
	//solve for dInvLHS
	GetInvLHS(B);
	//solve for dMu
	//EvalDensity(B);
	
	//make high res contour
	B->iHighResMult = 6;
	SZ = (size_t)(B->iHighResMult * (2*iTotalSize+3*iNLets*(NC-1)) ) * sizeof(double);
	B->dMuHigh = (double*)malloc(SZ);
	SZ = (size_t)NC * sizeof(Contour);
	B->CHigh = (Contour*) malloc(SZ);
	for (i=0; i<NC; i++)
	{
		N = B->iHighResMult * B->NP[i];
		Make_Contour (N, &(B->CHigh[i]) );
	}
	
	//upsample boundary, and curveprops
	UpsampleBoundary(B);
}

//calculates the inverse D matrix
void GetInvLHS(Boundary *B)
{
	int i, j, ind, N;
	int NC = B->NC; //number of contours
	int NP = B->iTotalSize; //total number of points
	int NI = NC - 1; //number of internal domains
	int *iPIndex;
	double *dD, *dDD;
	size_t SZ;
	double *dNX, *dNY; //outward normals
	double *dTX, *dTY; //Tangents
	double *dX, *dY; //X/Y points
	double *dSa, *dH, dHL; //Arclength and something
	double *dKap; //curvature
	int *iNDN; //some county thing for the outer contours
	int *iIndShuff; //another counter
	double *dRX, *dRY, *dRho; //length differences
	double *dCo; //some sort of coefficient
	double dTXX, dTXY, dTYX, dTYY; //tangent dot product
	double dSHK; //some coefficient
	double *dNxX, *dNxY; //scaled normals
	double *dNDNX, *dNDNY; //vectors that are sort of related to the normals
	int iNLets = B->iNLets; //I'm not sure what this is for
	int *iPiv;
	double *dLHS;
	double *dS, *dV; //the side and bottom of the LHS matrix
	int b, k, w;
	int iSC[3]; //indexes of dS
	int *iSR; //indexes of dV
	int iHead;
	double *dRow1, *dRow2; //two rows used for building dV
	double *dRowV;
	
	//make room for D
	SZ = (size_t)((2*NP)*(2*NP)) * sizeof(double);
	dD  = (double*)malloc(SZ);
	dDD = (double*)malloc(SZ);
	
	//make room for the other quantities
	SZ = (size_t)NP * sizeof(double);
	dNX = (double*) malloc(SZ);
	dNY = (double*) malloc(SZ);
	dTX = (double*) malloc(SZ);
	dTY = (double*) malloc(SZ);
	dX = (double*) malloc(SZ);
	dY = (double*) malloc(SZ);
	dSa = (double*) malloc(SZ);
	dH = (double*) malloc(SZ);
	dKap = (double*) malloc(SZ);
	
	//initialize iPIndex - starting positions of different contours
	SZ = (size_t)NC * sizeof(int);
	iPIndex = (int*)malloc(SZ);
	
	//copy in quantities from the contours
	ind = 0;
	for (i=0; i<NC; i++)
	{
		N = B->C[i].N;
		SZ = (size_t)N * sizeof(double);
		memcpy( &(dNX[ind]), B->C[i].dNormX, SZ);
		memcpy( &(dNY[ind]), B->C[i].dNormY, SZ);
		memcpy( &(dTX[ind]), B->C[i].dTangX, SZ);
		memcpy( &(dTY[ind]), B->C[i].dTangY, SZ);
		memcpy( &(dX[ind]), B->C[i].dXPts, SZ);
		memcpy( &(dY[ind]), B->C[i].dYPts, SZ);
		memcpy( &(dSa[ind]), B->C[i].dSa, SZ);
		memcpy( &(dKap[ind]), B->C[i].dKap, SZ);
		
		//h is the arc-length parameterization step size
		dHL = ( 2*M_PI )/( (double)N );
		for (j=0; j<N; j++)
			dH[ind+j] = dHL;

		ind += N;
		iPIndex[i] = ind;
	}
	
	//Initialize the counters
	SZ = (size_t)(2*B->C[0].N) * sizeof(int);
	iNDN = (int*) malloc(SZ);
	SZ = (size_t)(2*NP) * sizeof(int);
	iIndShuff = (int*) malloc(SZ);

	//fill in these funny counters
	N = B->C[0].N;
	for (i=0; i<N; i++)
	{
		iNDN[2*i] = i;
		iNDN[2*i + 1] = N+i;
	}
	for (i=0; i<NP; i++)
	{
		iIndShuff[i] = 2*i;
		iIndShuff[i+NP] = 2*i + 1;
	}
	
	//initialize distance counters, and the coefficients
	SZ = (size_t)NP * sizeof(double);
	dRX = (double*) malloc(SZ);
	dRY = (double*) malloc(SZ);
	dRho = (double*) malloc(SZ);
	dCo = (double*) malloc(SZ);

	//initialize these scaled normals
	SZ = (size_t)NP * sizeof(double);
	dNxX = (double*) malloc(SZ);
	dNxY = (double*) malloc(SZ);
	
	
	//initialize these guys
	SZ = (size_t)(2*NP) * sizeof(double);
	dNDNX = (double*) malloc(SZ);
	dNDNY = (double*) malloc(SZ);
	
	//double layer loop
	N = B->C[0].N; //number of points in outer contour
	for (i=0; i<NP; i++)
	{
		//get distances
		for (j=0; j<NP; j++)
		{
			dRX[j]  = dX[i]-dX[j];
			dRY[j]  = dY[i]-dY[j];
			dRho[j] = pow(dRX[j],2) + pow(dRY[j],2);
		}
		dRho[i] = 1; //ith Rho is 1, not 0
		
		//get coefficient
		for (j=0; j<NP; j++)
			dCo[j] = dSa[j] * dH[j] * ( dRX[j]*dNX[j] + dRY[j]*dNY[j] ) / ( pow(dRho[j],2) * M_PI );
		
		//fill in dD
		ind = (2*i)*(2*NP);
		for (j=0; j<NP; j++)
		{
			dD[ind+iIndShuff[j]]    = dCo[j]*dRX[j]*dRX[j];
			dD[ind+iIndShuff[j+NP]] = dCo[j]*dRX[j]*dRY[j];
		}
		ind = (2*i + 1)*(2*NP);
		for (j=0; j<NP; j++)
		{
			dD[ind+iIndShuff[j]]    = dCo[j]*dRY[j]*dRX[j];
			dD[ind+iIndShuff[j+NP]] = dCo[j]*dRY[j]*dRY[j];
		}
		
		//fill in dDD
		ind = i*(2*NP);
		for (j=0; j<NP; j++)
		{
			dDD[ind+j]    = dCo[j]*dRX[j]*dRX[j];
			dDD[ind+j+NP] = dCo[j]*dRX[j]*dRY[j];
		}
		ind = (i+NP)*(2*NP);
		for (j=0; j<NP; j++)
		{
			dDD[ind+j]    = dCo[j]*dRY[j]*dRX[j];
			dDD[ind+j+NP] = dCo[j]*dRY[j]*dRY[j];
		}
		
		//calculate tangent dot products
		dTXX = dTX[i] * dTX[i];
		dTXY = dTX[i] * dTY[i];
		dTYX = dTY[i] * dTX[i];
		dTYY = dTY[i] * dTY[i];
		
		//fill in special cases for dDD
		dSHK = dSa[i]*dH[i]*dKap[i] / (2*M_PI);
		dDD[i*(2*NP)      + i   ] = -dSHK*dTXX;
		dDD[i*(2*NP)      + i+NP] = -dSHK*dTXY;
		dDD[(i+NP)*(2*NP) + i   ] = -dSHK*dTYX;
		dDD[(i+NP)*(2*NP) + i+NP] = -dSHK*dTYY;
		
		//fill in special cases for dD
		dD[(2*i)*(2*NP)     + (2*i)    ] = -dSHK*dTXX;
		dD[(2*i)*(2*NP)     + (2*i + 1)] = -dSHK*dTXY;
		dD[(2*i + 1)*(2*NP) + (2*i)    ] = -dSHK*dTYX;
		dD[(2*i + 1)*(2*NP) + (2*i + 1)] = -dSHK*dTYY;		
		
		//special case for outer contour
		if (i < N)
		{
			for (j=0; j<N; j++)
			{
				dNxX[j] = dSa[j]*dH[j]*dNX[i];
				dNxY[j] = dSa[j]*dH[j]*dNY[i];
				dNDNX[j]   = dNxX[j]*dNX[j];
				dNDNX[j+N] = dNxX[j]*dNY[j];
				dNDNY[j]   = dNxY[j]*dNX[j];
				dNDNY[j+N] = dNxY[j]*dNY[j];
			}
			
			for (j=0; j<(2*N); j++)
			{			
				dD[(2*i)*(2*NP)   + j] += dNDNX[iNDN[j]];
				dD[(2*i+1)*(2*NP) + j] += dNDNY[iNDN[j]];
			}
		}
		
	}
	
	//subtract 0.5*I from dD
	N = 2*NP;
	for (i=0; i<N; i++)
		dD[i*N + i] -= 0.5;
	
	//for subdomains
	if (NI > 0) //at least 1 subdomain
	{
		//initialize S and V
		dS = (double*)calloc( (2*NP)*(3*iNLets*NI), sizeof(double));
		dV = (double*)calloc( (3*iNLets*NI)*(2*NP), sizeof(double));
		
		//get contributions from stokeslets and rotlets
		for (b=1; b<NC; b++)
			for (k=0; k<iNLets; k++)
			{
				//column indexes for dS
				iSC[0] = 2*(b-1)*iNLets - 2 + 2*(k+1);
				iSC[1] = 2*(b-1)*iNLets - 1 + 2*(k+1);
				
				//calculate distance from centerpoint
				for (i=0; i<NP; i++)
				{
					dRX[i]  = dX[i] - B->C[b].dCX;
					dRY[i]  = dY[i] - B->C[b].dCY;
					dRho[i] = pow(dRX[i], 2) + pow(dRY[i], 2);
				}
				
				//fill in S
				w = 3*iNLets*NI; //width of rows
				for (i=0; i<NP; i++)
				{
					dS[iIndShuff[i   ]*w + iSC[0]] = dRX[i]*dRX[i] / (4*M_PI*dRho[i]);
					dS[iIndShuff[i   ]*w + iSC[1]] = dRX[i]*dRY[i] / (4*M_PI*dRho[i]);
					dS[iIndShuff[i+NP]*w + iSC[0]] = dRY[i]*dRX[i] / (4*M_PI*dRho[i]);
					dS[iIndShuff[i+NP]*w + iSC[1]] = dRY[i]*dRY[i] / (4*M_PI*dRho[i]);
				}
				for (i=0; i<NP; i++)
				{
					dS[(2*i*w)     + iSC[0]] -= log(dRho[i])/(8*M_PI);
					dS[(2*i + 1)*w + iSC[1]] -= log(dRho[i])/(8*M_PI);
				}
				
				//sub column index for R
				iSC[2] = (b-1)*iNLets + k + 2*iNLets*NI;
				for (i=0; i<NP; i++)
				{
					dS[(2*i*w)     + iSC[2]] =  dRY[i] / dRho[i];
					dS[(2*i + 1)*w + iSC[2]] = -dRX[i] / dRho[i];
				}				
			}
			
		//make the right length for iSR
		SZ = (size_t)(2*iNLets) * sizeof(int);
		iSR = (int*)malloc(SZ);

		//prepare dRows
		N = B->C[1].N;
		dRow1 = (double*) calloc(2*N, sizeof(double));
		dRow2 = (double*) calloc(2*N, sizeof(double));
		SZ = (size_t)(2*N) * sizeof(double);
		dRowV = malloc(SZ);
		
		//calculate alpha and beta
		for (b=1; b<NC; b++)
		{
			N = B->C[b].N;
			
			//The index guy gets updated again
			for(i=0; i<2*iNLets; i++)
				iSR[i] = 2*iNLets*b-2*iNLets + i;
			
			//calculate head (pressure rise across the turbo-pump)
			iHead = 2*iPIndex[b-1];
			
			//realloc dRows
			SZ = (size_t)(2*N) * sizeof(double);
			dRow1 = (double*) realloc(dRow1, SZ);
			dRow2 = (double*) realloc(dRow2, SZ);
			dRowV = (double*) realloc(dRowV, SZ);
			
			//wipe dRow1 and dRow2
			memset(dRow1, 0, 2*N);
			memset(dRow2, 0, 2*N);
			
			//fill in row 1
			for (i=0; i<N; i++)
				dRow1[2*i] = dH[iPIndex[b-1]+i]*dSa[iPIndex[b-1]+i] / (2*M_PI*iNLets);
			
			//copy with offset to row 2
			SZ = (size_t)(2*N - 1) * sizeof(double);
			memcpy( &(dRow2[1]), &(dRow1[0]), SZ);
			
			//copy into dV
			w = (2*NP);
			SZ = (size_t)2*N * sizeof(double);
			for(i=0; i<iNLets; i++)
			{
				memcpy( &(dV[iSR[2*i    ]*w + iHead]), dRow1, SZ);
				memcpy( &(dV[iSR[2*i + 1]*w + iHead]), dRow2, SZ);
			}
			
			//get the next SR indices
			for (i=0; i<iNLets; i++)
				iSR[i] = 2*NI*iNLets + (b-1)*iNLets + i; 
			
			//refill dRow1
			for (i=0; i<N; i++)
				dRow1[i] = dH[iPIndex[b-1]+i]*dSa[iPIndex[b-1]+i] / (2*M_PI*iNLets);
			
			//now scale dRows by X,Y
			for (i=0; i<N; i++)
			{
				dRowV[2*i    ] = dRow1[i] *  B->C[b].dYPts[i];
				dRowV[2*i + 1] = dRow1[i] * -B->C[b].dXPts[i];
			}

			//copy into dV
			for(i=0; i<iNLets; i++)
				memcpy( &(dV[iSR[i]*w + iHead]), dRowV, SZ);
		}	
		
		//assemble the fat matrix
		N = 2*NP+3*iNLets*NI;
		SZ = (size_t)(N*N) * sizeof(double);
		dLHS = (double*) malloc(SZ);
		//copy in dD
		SZ = (size_t)(2*NP) * sizeof(double);
		for (i=0; i<(2*NP); i++)
			memcpy( &(dLHS[i*N]), &(dD[i*2*NP]), SZ);
		//copy in dS
		SZ = (size_t)(3*NI*iNLets) * sizeof(double);
		for (i=0; i<(2*NP); i++)
			memcpy( &(dLHS[i*N + 2*NP]), &(dS[i*3*NI*iNLets]), SZ);
		//copy in dV
		SZ = (size_t)(2*NP) * sizeof(double);
		for (i=0; i<(3*NI*iNLets); i++)
			memcpy( &(dLHS[(2*NP + i)*N]), &(dV[i*2*NP]), SZ);
		//stick -I in the bottom right corner
		for (i=0; i<(3*NI*iNLets); i++)
			for (j=0; j<(3*NI*iNLets); j++)
				dLHS[(2*NP+i)*N + 2*NP + j] = (i==j) ? -1 : 0;
		//local clean up
		free(dS); free(dV);
		free(iSR);
	}
	else //no subdomains
	{
		//copy dD into dLHS
		N = 2*NP;
		SZ = (size_t)(N*N) * sizeof(double);
		dLHS = (double*) malloc(SZ);
		memcpy(dLHS, dD, SZ);
	}
	
	//get matrix inverse
	SZ = (size_t)N * sizeof(int);
	iPiv = (int*) malloc(SZ);
	clapack_dgetrf(CblasRowMajor, N, N, dLHS, N, iPiv);
	clapack_dgetri(CblasRowMajor, N, dLHS, N, iPiv);
	
	//copy output
	SZ = (size_t)(N*N) * sizeof(double);
	memcpy(B->dInvLHS, dLHS, SZ);
	
	//clean up
	free(iPIndex);
	free(dD); free(dDD);
	free(dNX); free(dNY);
	free(dTX); free(dTY);
	free(dX); free(dY);
	free(dSa); free(dH);
	free(dKap);
	free(iNDN); free(iIndShuff);
	free(dRX); free(dRY); free(dRho);
	free(dCo);
	free(dNxX); free(dNxY);
	free(dNDNX); free(dNDNY);
	free(iPiv); free(dLHS);
	free(dRow1); free(dRow2);
	free(dRowV);
}

//uses upsampled boundary to calculate dU/dV at dX/dY
//This is only good for iNLets=1
void EvalBDVelocity(int NIn, double *dXIn, double *dYIn, double *dU, double *dV, int iMode, ...)
{
	size_t SZ;
	int NP; 
	int NC;
	int NI; 
	int iNLets;
	int i, j, b, k;
	int ind, N;
	//joined curve props
	double *dX, *dY; //X/Y points
	double *dNX, *dNY; //outward normals
	double *dTX, *dTY; //Tangents
	double *dSa, *dH, dHL; //Arclength and arc length step size
	double *dKap; //curvature
	double *dCX, *dCY; //centerpoints
	double *dA; //matrix for solving velocities
	int *iIndShuff;
	double *dRX, *dRY, *dRho;
	double *dCo;
	double *dS, *dR; //stokeslet and rotlet contributions
	double *dSH; //combines Sa and H to save time
	int iSC[3];
	int m, n;
	double *dYou; //output vector
	static Boundary *B;
	va_list vaArgs;
	
	if (iMode == 1) //initialization
	{
		va_start (vaArgs, iMode);    
		B = va_arg(vaArgs, Boundary*);
		va_end(vaArgs);
	}
	else //regular velocity calculation
	{
		//set these numbers
		NP = B->iTotalSize * B->iHighResMult;
		NC = B->NC;
		NI = NC - 1;
		iNLets = B->iNLets;
		
		//initialize big curve props vectors
		SZ = (size_t)NP * sizeof(double);
		dX = (double*)malloc(SZ);
		dY = (double*)malloc(SZ);
		dNX = (double*)malloc(SZ);
		dNY = (double*)malloc(SZ);
		dTX = (double*)malloc(SZ);
		dTY = (double*)malloc(SZ);
		dSa = (double*)malloc(SZ);
		dH = (double*)malloc(SZ);
		dKap = (double*)malloc(SZ);
		dSH = (double*)malloc(SZ);
		SZ = (size_t)NC * sizeof(double);
		dCX = (double*)malloc(SZ);
		dCY = (double*)malloc(SZ);
		
		//join up curve props
		ind = 0;
		for (i=0; i<NC; i++)
		{
			N = B->CHigh[i].N;
			SZ = (size_t)N * sizeof(double);
			
			memcpy( &(dNX[ind]),  B->CHigh[i].dNormX, SZ);
			memcpy( &(dNY[ind]),  B->CHigh[i].dNormY, SZ);
			memcpy( &(dTX[ind]),  B->CHigh[i].dTangX, SZ);
			memcpy( &(dTY[ind]),  B->CHigh[i].dTangY, SZ);
			memcpy( &(dX[ind]),   B->CHigh[i].dXPts, SZ);
			memcpy( &(dY[ind]),   B->CHigh[i].dYPts, SZ);
			memcpy( &(dSa[ind]),  B->CHigh[i].dSa, SZ);
			memcpy( &(dKap[ind]), B->CHigh[i].dKap, SZ);
			
			//h is the arc-length parameterization step size
			dHL = ( 2*M_PI )/( (double)N );
			for (j=0; j<N; j++)
			{
				dH[ind+j] = dHL;
				dSH[ind+j] = dHL*dSa[ind+j]/(M_PI);
			}
			
			dCX[i] = B->CHigh[i].dCX;
			dCY[i] = B->CHigh[i].dCY;
			
			ind += N;
		}
		
		//initialize main matrix, and iIndShuff
		SZ = (size_t)(2*NIn  * 2*NP) * sizeof(double);
		dA = (double*)malloc(SZ);
		SZ = (size_t)(2*NP) * sizeof(int);
		iIndShuff = (int*)malloc(SZ);
		
		//set iIndShuff
		for (i=0; i<NP; i++)
		{
			iIndShuff[i   ] = 2*i;
			iIndShuff[i+NP] = 2*i + 1;	
		}
		
		//initialize distance measures
		SZ = (size_t)NP * sizeof(double);
		dRX = (double*)malloc(SZ);
		dRY = (double*)malloc(SZ);
		dRho = (double*)malloc(SZ);
		dCo = (double*)malloc(SZ);
		
		//fill in main block of A
		for (i=0; i<NIn; i++)
		{
			//calculate distances
			for (j=0; j<NP; j++)
			{
				dRX[j] = dXIn[i] - dX[j];
				dRY[j] = dYIn[i] - dY[j];
				dRho[j] = pow(dRX[j],2) + pow(dRY[j],2);
				dCo[j] = dSH[j] * (dRX[j]*dNX[j] + dRY[j]*dNY[j]) / ( dRho[j]*dRho[j] );
			}
			
			//fill in entries of A
			for (j=0; j<NP; j++)
			{
				dA[(2*i)*(2*NP)     + iIndShuff[j]   ] = dCo[j]*dRX[j]*dRX[j];
				dA[(2*i)*(2*NP)     + iIndShuff[j+NP]] = dCo[j]*dRX[j]*dRY[j];
				dA[(2*i + 1)*(2*NP) + iIndShuff[j]   ] = dCo[j]*dRY[j]*dRX[j];
				dA[(2*i + 1)*(2*NP) + iIndShuff[j+NP]] = dCo[j]*dRY[j]*dRY[j];
			}
		}
		
		//initialize dS, dR
		dS = (double*) calloc(2*NIn * 2*NI*iNLets, sizeof(double));
		dR = (double*) calloc(2*NIn * NI*iNLets, sizeof(double));
		
		//resize dRX, dRY, dRho
		SZ = (size_t)(NIn) * sizeof(double);
		dRX = (double*)realloc(dRX, SZ);
		dRY = (double*)realloc(dRY, SZ);
		dRho = (double*)realloc(dRho, SZ);
		
		//add stokeslet/rotlet contributions
		for (b=1; b<NC; b++)
			for (k=0; k < iNLets; k++)
			{
				//column indexes for dS
				iSC[0] = 2*(b-1)*iNLets - 2 + 2*(k+1);
				iSC[1] = 2*(b-1)*iNLets - 1 + 2*(k+1);
				
				//fill in dRX/dRY/dRho
				for (i=0; i<NIn; i++)
				{
					//this only supports single stokeslet/rotlets
					dRX[i] = dXIn[i] - dCX[b];
					dRY[i] = dYIn[i] - dCY[b];
					dRho[i] = pow(dRX[i],2) + pow(dRY[i],2);
				}
				
				//fill in dS
				for (i=0; i<NIn; i++)
				{
					dS[(2*i)    *(2*NI*iNLets) + iSC[0]] = dRX[i]*dRX[i]/(4*M_PI*dRho[i]);
					dS[(2*i + 1)*(2*NI*iNLets) + iSC[0]] = dRX[i]*dRY[i]/(4*M_PI*dRho[i]);
					dS[(2*i)    *(2*NI*iNLets) + iSC[1]] = dRY[i]*dRX[i]/(4*M_PI*dRho[i]);
					dS[(2*i + 1)*(2*NI*iNLets) + iSC[1]] = dRY[i]*dRY[i]/(4*M_PI*dRho[i]);
					dS[(2*i)    *(2*NI*iNLets) + iSC[0]] -= log(dRho[i]) / (8 * M_PI);
					dS[(2*i + 1)*(2*NI*iNLets) + iSC[1]] -= log(dRho[i]) / (8 * M_PI);
				}
			
				//fill in dR
				iSC[2] = (b-1)*iNLets + k;
				for (i=0; i<NIn; i++)
				{
					dR[(2*i    )*NI*iNLets + iSC[2]] =  dRY[i]/dRho[i];
					dR[(2*i + 1)*NI*iNLets + iSC[2]] = -dRX[i]/dRho[i];
				}	
			}
		
		
		//get velocity vector, contribution from A
		SZ = (size_t)(2*NIn) * sizeof(double);
		dYou = (double*)malloc(SZ);
		m = 2*NIn;
		n = 2*NP;
		cblas_dgemv(CblasRowMajor, CblasNoTrans, m, n, 1.0, dA, n, B->dMuHigh, 1, 0.0, dYou, 1);
		
		if (NI > 0) //add contribution from S and R
		{
			n = 2*NI;
			i = 2*NP;
			cblas_dgemv(CblasRowMajor, CblasNoTrans, m, n, 1.0, dS, n, &(B->dMuHigh[i]), 1, 1, dYou, 1);
			n = NI;
			i += 2*NI;
			cblas_dgemv(CblasRowMajor, CblasNoTrans, m, n, 1.0, dR, n, &(B->dMuHigh[i]), 1, 1, dYou, 1);
		}
		
		//stick it in dU/dV
		for (i=0; i<NIn; i++)
		{
			dU[i] = dYou[2*i];
			dV[i] = dYou[2*i + 1];
		}
		
		//clean up
		free(dX); free(dY);
		free(dNX); free(dNY);
		free(dTX); free(dTY);
		free(dSa); free(dH);
		free (dKap);
		free(dCX); free(dCY);
		free(dA);
		free(iIndShuff);
		free(dRX); free(dRY); free(dRho);
		free(dCo);
		free(dS); free(dR); 
		free(dYou);
		free(dSH);
	}
}

void UpsampleBoundary(Boundary *B)
{
	int NC = B->NC;
	int N, NH;
	int i;
		
	//upsample all of the curve props
	for (i=0; i<NC; i++)
	{
		N = B->NP[i];
		NH = N *B->iHighResMult;
		B->CHigh[i].N = NH;
		Interp_Double(B->C[i].dXPts, N, B->CHigh[i].dXPts, NH);
		Interp_Double(B->C[i].dYPts, N, B->CHigh[i].dYPts, NH);
		Interp_Double(B->C[i].dU, N, B->CHigh[i].dU, NH);
		Interp_Double(B->C[i].dV, N, B->CHigh[i].dV, NH);
		Interp_Double(B->C[i].dDX, N, B->CHigh[i].dDX, NH);
		Interp_Double(B->C[i].dDY, N, B->CHigh[i].dDY, NH);
		Interp_Double(B->C[i].dDXISa, N, B->CHigh[i].dDXISa, NH);
		Interp_Double(B->C[i].dDYISa, N, B->CHigh[i].dDYISa, NH);
		Interp_Double(B->C[i].dDDX, N, B->CHigh[i].dDDX, NH);
		Interp_Double(B->C[i].dDDY, N, B->CHigh[i].dDDY, NH);
		Interp_Double(B->C[i].dD4XISa, N, B->CHigh[i].dD4XISa, NH);
		Interp_Double(B->C[i].dD4YISa, N, B->CHigh[i].dD4YISa, NH);
		Interp_Double(B->C[i].dTangX, N, B->CHigh[i].dTangX, NH);
		Interp_Double(B->C[i].dTangY, N, B->CHigh[i].dTangY, NH);
		Interp_Double(B->C[i].dNormX, N, B->CHigh[i].dNormX, NH);
		Interp_Double(B->C[i].dNormY, N, B->CHigh[i].dNormY, NH);
		Interp_Double(B->C[i].dKap, N, B->CHigh[i].dKap, NH);
		Interp_Double(B->C[i].dSa, N, B->CHigh[i].dSa, NH);
		Interp_Double(B->C[i].dISa, N, B->CHigh[i].dISa, NH);
		B->CHigh[i].dCX = B->C[i].dCX;
		B->CHigh[i].dCY = B->C[i].dCY;
	}
}

void UpsampleMu(Boundary *B)
{
	int NC = B->NC;
	int i, j, ind, indHigh, N, NH;
	size_t SZ;
	double *dMuX, *dMuY;
	double *dMuXHigh, *dMuYHigh;

	//Initialize dMuX/Y and the high res versions
	SZ = (size_t)B->NP[0] * sizeof(double);
	dMuX = (double*)malloc(SZ);
	dMuY = (double*)malloc(SZ);
	SZ *= (size_t)B->iHighResMult;
	dMuXHigh = (double*)malloc(SZ);
	dMuYHigh = (double*)malloc(SZ);
	
	//upsample dMu
	ind = 0; indHigh = 0;
	for (i=0; i<NC; i++)
	{
		N = B->NP[i];
		NH = N *B->iHighResMult;
		
		//make sure dMus are the right size
		SZ = (size_t)N * sizeof(double);
		dMuX = (double*)realloc(dMuX, SZ);
		dMuY = (double*)realloc(dMuY, SZ);
		SZ *= (size_t)NH * sizeof(double);
		dMuXHigh = (double*)realloc(dMuXHigh, SZ);
		dMuYHigh = (double*)realloc(dMuYHigh, SZ);
		
		//shuffle in the old dMu values
		for (j=0; j<N; j++)
		{
			dMuX[j] = B->dMu[2*j + ind];
			dMuY[j] = B->dMu[2*j + ind + 1];
		}
		
		//upsample them into dMuX/YHigh
		Interp_Double(dMuX, N, dMuXHigh, NH);
		Interp_Double(dMuY, N, dMuYHigh, NH);
		
		//shuffle high res into dMuHigh
		for (j=0; j<NH; j++)
		{
			B->dMuHigh[indHigh + 2*j]     = dMuXHigh[j];
			B->dMuHigh[indHigh + 2*j + 1] = dMuYHigh[j];
		}
		
		//update indices
		ind += (2*N);
		indHigh += (2*NH);
	}
	
	//copy in the stokeslet/rotlet contributions
	SZ = (size_t)( 3*B->iNLets*(NC-1) ) * sizeof(double);
	memcpy( &(B->dMuHigh[indHigh]), &(B->dMu[ind]), SZ);
	
	//clean up
	free(dMuX); free(dMuY);
	free(dMuXHigh); free(dMuYHigh);
}

void Dest_Boundary(Boundary *B)
{
	int i;
	
	for(i=0; i<(B->NC); i++)
	{	
		Dest_Contour( &(B->C[i]) );
		Dest_Contour( &(B->CHigh[i]) );
	}
	free(B->CHigh);
	free(B->dInvLHS);
	free(B->NP);
	free(B->dMu);
	free(B->dMuHigh);
	free(B->C);
}
