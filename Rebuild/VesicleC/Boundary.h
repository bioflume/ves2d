/*******************************************************************************
	Boundary.h															7/24/09
	Alexander S Rattner

Handles all of the fixed boundary specific code.
*******************************************************************************/

#ifndef BOUNDARY_H
#define BOUNDARY_H

//bring in externals
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cblas.h>
#include <clapack.h>
#include <stdarg.h>
#include <math.h>
#include "Contour.h"
#include "Misc.h"

//This contains all of the information that describes a boundary
typedef struct
{
	int NC; //the number of contours in the boundary
	int *NP; //this is a convenient array of the number of points in each contour
	int iTotalSize; //The total number of points
	//The outer contour should be ccw, the inner ones should be cw
	Contour *C; //the array of contours
	double *dInvLHS; //This matrix is the operator that solves for velocity inside the boundary
	int iNLets; //Number of Stokeslet/Rotlets in each subdomain
	double *dMu; //the density vector
	int iHighResMult; //the high resolution multiplier (NPHigh = iHighResMultiplier*NP)
	Contour *CHigh; //high resolution contour
	double *dMuHigh; //high resolution density
} Boundary;

void Make_Boundary(int iNLets, int NC, Contour *C, Boundary *B);
void GetInvLHS(Boundary *B);
void EvalBDVelocity(int NIn, double *dXIn, double *dYIn, double *dU, double *dV, int iMode, ...);
void UpsampleBoundary(Boundary *B);
void UpsampleMu(Boundary *B);
void Dest_Boundary(Boundary *B);

#endif
