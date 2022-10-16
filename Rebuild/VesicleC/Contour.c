/*******************************************************************************
	Contour.h														8/17/09
	Alexander S Rattner

Defines a closed 2D contour and helper functions
*******************************************************************************/

#include "Contour.h"


//Maker functions
void Make_Contour (int N, Contour *C)
{
	size_t SZ;
	C->N = N;
	
	//Allocate the memory for the boundary (and its velocities) and assign it to the vesicle
	SZ = (size_t)N*sizeof(double);
	C->dXPts = (double*) malloc(SZ); //Puts the pointer to dXpts
	C->dYPts = (double*) malloc(SZ);
	//I suppose that you may want to give the vesicles an initial velocity
	C->dU    = (double*) calloc(N, sizeof(double));
	C->dV    = (double*) calloc(N, sizeof(double));
	
	//Malloc for derivatives
	C->dDX = (double*) malloc(SZ);
	C->dDY = (double*) malloc(SZ);
	C->dDXISa = (double*) malloc(SZ);
	C->dDYISa = (double*) malloc(SZ);
	C->dDDX = (double*) malloc(SZ);
	C->dDDY = (double*) malloc(SZ);
	C->dD4XISa = (double*) malloc(SZ);
	C->dD4YISa = (double*) malloc(SZ);
	
	//Malloc for curve properties
	C->dTangX = (double*) malloc(SZ);
	C->dTangY = (double*) malloc(SZ);
	C->dNormX = (double*) malloc(SZ);
	C->dNormY = (double*) malloc(SZ);
	C->dKap = (double*) malloc(SZ);
	C->dSa = (double*) malloc(SZ);
	C->dISa = (double*) malloc(SZ);
}

//Initializes a simple contour data structure
SimContour Make_SimContour(int N)
{
	SimContour SC;
	size_t SZ;
	
	SC.N = N;
	
	SZ = (size_t)N * sizeof(double);
	SC.dXPts = (double*) malloc(SZ);
	SC.dYPts = (double*) malloc(SZ);
	
	return SC;
}

//Gets the curve properties from a Contour
//I assume that the calling function malloced the curve props for me
void CurveProps(Contour *C)
{
	int i;
	int N = C->N;
	double dDist;
	
	//Compute regular derivatives
	DmFT_Double (N, 1, C->dXPts, C->dDX, 0);
	DmFT_Double (N, 1, C->dYPts, C->dDY, 0);
	DmFT_Double (N, 1, C->dDX, C->dDDX, 0);
	DmFT_Double (N, 1, C->dDY, C->dDDY, 0);
	
	//Compute curve props
	for (i=0; i<N; i++)
	{
		//calculate arc length (dSa)
		C->dSa[i] = sqrt( C->dDX[i]*C->dDX[i] + C->dDY[i]*C->dDY[i] );
		//This is a holdover from the MATLAB code, it looks pretty hackish
		if (C->dSa[i] == 0)
			C->dSa[i] = 1;
		
		//get inverse arc length
		C->dISa[i] = 1 / C->dSa[i];

		//calculate local tangent (dTangX, dTangY)
		C->dTangX[i] = C->dDX[i]/C->dSa[i];
		C->dTangY[i] = C->dDY[i]/C->dSa[i];
		
		//calculate local normals (outward)
		C->dNormX[i] =  C->dTangY[i];
		C->dNormY[i] = -C->dTangX[i];
		
		//calculate kurvature (dKap)
		C->dKap[i] = (C->dDX[i]*C->dDDY[i] - C->dDY[i]*C->dDDX[i]) / pow(C->dSa[i], 3);
	}
	
	//Compute scaled derivatives
	DmFT_Double (N, 1, C->dXPts, C->dDXISa, 1, C->dISa);
	DmFT_Double (N, 1, C->dYPts, C->dDYISa, 1, C->dISa);
	DmFT_Double (N, 3, C->dDXISa, C->dD4XISa, 1, C->dISa);
	DmFT_Double (N, 3, C->dDYISa, C->dD4YISa, 1, C->dISa);
	
	//calculate center point
	//this isn't necessarily the center point
	C->dCX = 0; C->dCY = 0;
	for (i=0; i<N; i++)
	{
		C->dCX += C->dXPts[i];
		C->dCY += C->dYPts[i];
	}
	C->dCX /= (double)N;
	C->dCY /= (double)N;
	
	//calculate maximum distance between any 2 points
	C->dMaxDist = 0;
	for (i=0; i<N; i++)
	{
		//get distance
		if (i<N-1)
		{	dDist = sqrt( pow(C->dXPts[i]- C->dXPts[i+1], 2) + pow(C->dYPts[i]- C->dYPts[i+1], 2) ); }
		else //last point
		{	dDist = sqrt( pow(C->dXPts[i]- C->dXPts[0], 2) + pow(C->dYPts[i]- C->dYPts[0], 2) ); }
		
		//check for maximum distance
		if (dDist > C->dMaxDist)
			C->dMaxDist = dDist;
	}
}

//Clears Vesicle from Memory
void Dest_Contour(Contour *C)
{
	free(C->dXPts); free(C->dYPts);
	free(C->dU); free(C->dV);
	free(C->dDX); free(C->dDY);
	free(C->dDXISa); free(C->dDYISa);
	free(C->dDDX); free(C->dDDY);
	free(C->dD4XISa); free(C->dD4YISa);
	free(C->dTangX); free(C->dTangY);
	free(C->dNormX); free(C->dNormY);
	free(C->dKap); 
	free(C->dSa); free(C->dISa);
}

void Dest_SimContour(SimContour SC)
{
	free(SC.dXPts);
	free(SC.dYPts);
}
