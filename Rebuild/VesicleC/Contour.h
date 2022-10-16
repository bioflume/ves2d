/*******************************************************************************
	Contour.h														8/17/09
	Alexander S Rattner

Defines a closed 2D contour and helper functions
*******************************************************************************/

#ifndef CONTOUR_H
#define CONTOUR_H

//Pull in externals
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <cblas.h>
#include <fftw3.h>
#include "Misc.h"

//Define Constants here
#define MAX_DESC_LENGTH 30

typedef struct
{
	int N; //The number of points on the vesicle boundary
	double *dXPts, *dYPts; //X, Y points on vesicle boundary
	double *dU, *dV; //X,Y velocities on boundary
	double *dDX, *dDY; //First derivatives
	double *dDXISa, *dDYISa; //First derivative, scaled by inverse arc length
	double *dDDX, *dDDY; //Second derivatives
	double *dD4XISa, *dD4YISa; //4th derivatives
	double *dTangX, *dTangY; //curve tangents
	double *dNormX, *dNormY;
	double *dKap; //curvature
	double *dSa, *dISa; //arclength and inverse arc-length
	double dCX, dCY; //center points
	double dMaxDist; //maximum distance between two points
} Contour;

//A simplified contour
typedef struct
{
	int N;
	double *dXPts, *dYPts;
} SimContour;

//Makers
void Make_Contour (int N, Contour *C);
SimContour Make_SimContour(int N);

//other functions
void CurveProps(Contour *C);

//Destroyers
void Dest_Contour(Contour *C);
void Dest_SimContour(SimContour SC);

#endif
