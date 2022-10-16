/*******************************************************************************
	Misc.c															6/25/09
	Alexander S Rattner

Does some auxillary stuff
*******************************************************************************/

#include "Misc.h"

//performs an fftshift on array F of size SZ
void fftshift(void *F, size_t SZ)
{
	void *Temp;
	Temp = malloc(SZ/2);
	
	//swap first and second halves
	memcpy(Temp, F, (SZ/2));
	memcpy(F,F+(SZ/2), (SZ/2));
	memcpy(F+(SZ/2), Temp, (SZ/2));
	
	free(Temp);
}

//takes the mth derivative of an N length real function
//pF is the forward transform, pB is the backwards transform
//cF is the output of the forward transform, the derivative is stored in the pB output
void DmFT_R(int N, int m, double *dIn, fftw_plan pF, fftw_complex *cF, fftw_plan pB, double *dOut, int iMode, ...)
{
	int i, j, l;
	size_t SZ;
	fftw_complex cI, cIm;
	double *dS;
	double *dInStore;
	fftw_complex *cK;
	va_list vaArgs;
	
	//back up input vector, for higher orders
	SZ = (size_t)N * sizeof(double);
	if (m > 1)
	{
		dInStore = (double*) malloc(SZ);
		memcpy(dInStore, dIn, SZ);
	}
	
	//l is the actual output length (not the logical length)
	l = (N/2)+1;
	
	//prepare I^m
	cI[0] = 0; cI[1] = 1;
	CompPow(cI, m, cIm);
	
	//handle multiple arguments
	va_start (vaArgs, iMode);    
	
	//handle cK
	cK = (fftw_complex*) calloc(l, sizeof(fftw_complex));
	for (i=0; i<(l-1); i++)
		cK[i][1] = (double)i;
	cK[l-1][1] = -(N/2);

	//handle dS
	SZ = (size_t)N * sizeof(double);
	dS = (double*) malloc(SZ);
	if (iMode == 1)
	{	memcpy(dS, va_arg(vaArgs, double*), SZ); }
	else
	{
		for (i=0; i<N; i++)
			dS[i] = 1;
	}
	
	for (j=0; j<m; j++) //once for each order of derivative
	{
		//First compute the forward transform
		fftw_execute(pF);
		
		//Apply scaling factor
		for (i=0; i<l; i++)
			CompMult(cF[i], cK[i], cF[i]);
		
		//perform inverse transform
		fftw_execute(pB);
		
		//scale output
		for (i=0; i<N; i++)
		{
			dOut[i] = dS[i]*dOut[i] / (double)N;
		}
		
		//Stick it back in the input for higher orders
		if (j < (m-1))
			memcpy(dIn, dOut, SZ);
	}
	
	//restore original dIn for higher orders
	if (m > 1)
	{
		memcpy(dIn, dInStore, SZ);
		free(dInStore);
	}
	
	//clean up
	va_end(vaArgs);
	fftw_free(cK);
	free(dS);
}

//This is a simplified call to DmFT, it makes all of the fftw objects for you
void DmFT_Double (int N, int m, double *dIn, double *dOut, int iMode, ...)
{
	size_t SZ;
	double* dInP;
	double *dSa;
	fftw_plan pF, pB;
	fftw_complex *cF;
	
	
	//handle multiple arguments
	va_list vaArgs;
	va_start (vaArgs, iMode);    
	
	//make room for cF
	SZ = (size_t)N * sizeof(fftw_complex);
	cF = (fftw_complex*) fftw_malloc(SZ);
	
	//make room for dInP, we can't use dIn directly because the plan wipes it!
	SZ = (size_t)N * sizeof(double);
	dInP = (double*) malloc(SZ);
	
	//Make FFTW plans
	pF = fftw_plan_dft_r2c_1d(N, dInP, cF, FFTW_MEASURE | FFTW_PRESERVE_INPUT);
	pB = fftw_plan_dft_c2r_1d(N, cF, dOut, FFTW_MEASURE | FFTW_PRESERVE_INPUT);
	
	//copy dIn to dInP for above reason
	memcpy(dInP, dIn, SZ);
	
	switch(iMode)
	{
		case 0:
			DmFT_R(N, m, dInP, pF, cF, pB, dOut, 0);
			break;
		case 1:
			dSa = va_arg(vaArgs, double*);
			DmFT_R(N, m, dInP, pF, cF, pB, dOut, 1, dSa);
			break;
	}
	
	//clean up
	va_end(vaArgs);
	free(dInP);
	fftw_free(cF);
	fftw_destroy_plan(pF);
	fftw_destroy_plan(pB);
}

//does an fft interpolation of dX (length iX) to get dY (length (iY)
//dY should be longer than dX, and both should have even length
void Interp_Double(double *dX, int iX, double *dY, int iY)
{
	size_t SZ;
	fftw_complex *cX, *cY;	
	fftw_plan pF, pB;
	int i;
	int iNyquist;
	
	//initialize cX, cY, they should be iX, and iY long
	SZ = (size_t)iX * sizeof(fftw_complex);
	cX = (fftw_complex*) fftw_malloc(SZ);
	SZ = (size_t)iY * sizeof(fftw_complex);
	cY = (fftw_complex*) fftw_malloc(SZ);
	
	//make forward plans
	pF = fftw_plan_dft_1d(iX, cX, cY, FFTW_FORWARD, FFTW_MEASURE);
	pB = fftw_plan_dft_1d(iY, cY, cY, FFTW_BACKWARD, FFTW_MEASURE);

	//copy in iX to cX
	for (i=0; i<iX; i++)
	{	cX[i][0] = dX[i]; cX[i][1] = 0;}

	//run forward transform
	fftw_execute(pF);
	
	//reshape cY
	iNyquist = floor((iX+1)/2);
	SZ = (size_t)(iNyquist-1)*sizeof(fftw_complex);
	memcpy(&(cY[iY-iNyquist+1]), &(cY[iNyquist+1]), SZ);
	
	//clear out middle
	SZ = (size_t)(iY-iX) * sizeof(fftw_complex);
	memset( &(cY[iNyquist+1]), 0, SZ);

	//scale these two entries
	cY[iNyquist][0] /= 2; cY[iNyquist][1] /= 2;
	cY[iNyquist+iY-iX][0] = cY[iNyquist][0];
	cY[iNyquist+iY-iX][1] = cY[iNyquist][1];
	
	//run inverse transform
	fftw_execute(pB);
	
	//copy into dY
	for (i=0; i<iY; i++)
		dY[i] = cY[i][0];
		
	//scale output
	cblas_dscal(iY, 1 / (double)iX, dY, 1);
	
	//free up memory
	fftw_destroy_plan(pF);
	fftw_destroy_plan(pB);
	fftw_free(cX);
	fftw_free(cY);
}

//This creates the real interpolation matrix of dY in blas format
//This multiplies length iY dY to get length iX interpolant
double *GetInterpMat(double *dY, int iY, int iX, double dL)
{
	double dC, *dA;
	fftw_complex *cF, *cf, *cTemp; //FFT of dF
	fftw_plan p;
	int i, j, k;
	size_t szHalf, szFull;
	
	//dA is the output matrix
	dA = (double*) calloc(iX*iY, sizeof(double));
	
	//dC is some sort of scaling val
	dC = 2*M_PI/dL;
	
	//Make space for dF
	szFull = (size_t)iX*sizeof(fftw_complex);
	szHalf = (size_t)(iX/2)*sizeof(fftw_complex);
	cf = (fftw_complex*) fftw_malloc(szFull);
	cF = (fftw_complex*) fftw_malloc(szFull);
	cTemp = (fftw_complex*) fftw_malloc(szHalf);
	
	//clear cf
	for (i=0; i<iX; i++)
	{ cf[i][0] = 0; cf[i][1] = 0; }
	
	//prepare fft plan
	p = fftw_plan_dft_1d(iX, cf, cF, FFTW_FORWARD, FFTW_MEASURE);

	for (j=0; j<iX; j++)
	{
		//set f to zeros except at j
		if (j == 0)
		{	cf[j][0] = 1;}
		else
		{	cf[j-1][0] = 0; cf[j][0] = 1;}
		
		//run fft
		fftw_execute(p);
		//Scale output
		for (k=0; k<iX; k++)
		{
			cF[k][0] =  cF[k][0] / iX;
			cF[k][1] =  cF[k][1] / iX;
		}
		//Perform fftshift
		fftshift(cF, szFull);
		
		//This should be equivalent to the MATLAB expression (I hacked the Re(c_1*c_2) )
		for (k=0; k<(iX); k++)
			for (i=0; i<iY; i++)
				dA[(i*iX)+j] = dA[(i*iX)+j] + cF[k][0]*cos(dC*((-iX/2)+k)*dY[i])
										    - cF[k][1]*sin(dC*((-iX/2)+k)*dY[i]);	
		
	}

	//return memory
	fftw_free(cf); fftw_free(cF);
	fftw_destroy_plan(p);
	fftw_free(cTemp);
	
	return dA;
}

//Checks if test points are inside the polygon, adapted from W. Randolph Franklin's website
//(I couldn't find any copyright notice, so I assume it is ok to use this for general purposes)
//this treats boundary points as being inside (opposite of matlab)
void InPoly(int iNVert, double *dVertX, double *dVertY, int iNTest, double *dTestX, double *dTestY, int *iIn)
{
  int i, j, k;
	
	for (k=0; k<iNTest; k++)
	{
		iIn[k] = 0;
		for (i = 0, j = iNVert-1; i < iNVert; j = i++)
		{
			if ( ((dVertY[i]>dTestY[k]) != (dVertY[j]>dTestY[k])) &&
				 (dTestX[k] < (dVertX[j]-dVertX[i]) * (dTestY[k]-dVertY[i]) / (dVertY[j]-dVertY[i]) + dVertX[i]) )
				iIn[k] = !iIn[k];
		}
	}
}

//I don't know why this isn't built into math.h
double sign(double dIn)
{
	return ( (dIn >= 0) ? 1.0 : -1.0 );
}

//reads a string of numbers into an array of doubles
//if there are less than N numbers in str, this only fills in that many entries of dOut
void ParseArray(char *str, double *dOut, int N)
{
	int i = 0;
	int j = 0;
	int	nc;
	double dNum;
	
	while ( (j<N) && (sscanf( &(str[i]), "%lf%n", &dNum, &nc) == 1) )
	{
			dOut[j++] = dNum;
			i += nc;
	}
}

//simple complex power algorithm
inline void CompPow(fftw_complex cIn, int m, fftw_complex cOut)
{
	int i;	
	fftw_complex cTemp;
	
	cTemp[0] = cIn[0]; cTemp[1] = cIn[1];
	for (i=1; i<m; i++)
		CompMult(cIn, cTemp, cTemp);
	cOut[0] = cTemp[0]; cOut[1] = cTemp[1];
}

//simple complex multiplication algorithm
inline void CompMult(fftw_complex A, fftw_complex B, fftw_complex C)
{
	fftw_complex cTemp;
	cTemp[0] = A[0]*B[0] - A[1]*B[1];
	cTemp[1] = A[1]*B[0] + A[0]*B[1];
	C[0] = cTemp[0]; C[1] = cTemp[1];
}
