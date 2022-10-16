/*******************************************************************************
	Misc.h															6/25/09
	Alexander S Rattner

Does some auxillary stuff
*******************************************************************************/

#ifndef MISC_H
#define MISC_H

//Bring in externals
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>
#include <fftw3.h>
#include <stdarg.h>
#include <cblas.h>

//function prototypes
void fftshift(void *F, size_t SZ);
void DmFT_Double (int N, int m, double *dIn, double *dOut, int iMode, ...);
void DmFT_R(int N, int m, double *dIn, fftw_plan pF, fftw_complex *cF, fftw_plan pB, double *dOut, int iMode, ...);
void Interp_Double(double *dX, int iX, double *dY, int iY);
double *GetInterpMat(double *dY, int iY, int iX, double dL);
void InPoly(int iNVert, double *dVertX, double *dVertY, int iNTest, double *dTestX, double *dTestY, int *iIn);
double sign(double dIn);
void ParseArray(char *str, double *dOut, int N);

//Complex math stuff
inline void CompMult(fftw_complex A, fftw_complex B, fftw_complex C);
inline void CompPow(fftw_complex cIn, int m, fftw_complex cOut);

#endif
