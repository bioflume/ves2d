/*******************************************************************************
	StokesStructs.h														6/27/09
	Alexander S Rattner

Defines all of the structs needed in Stokes (Stokes is getting bloated) and 
their make/dest functions
*******************************************************************************/

#ifndef STOKESSTRUCTS_H
#define STOKESSTRUCTS_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include "f2c.h"
#include "Contour.h"
#include "Vesicle.h"
#include "DatOps.h"

//stores parameters used in GMRES
typedef struct
{
	integer iRestart;
	integer iLdw;
	integer iLdh;
	integer iIterMax;
	integer iIterAct;
	doublereal dResMax;
	doublereal dResAct;
	integer iInfo;
	int (*matvec) (doublereal*, doublereal*, doublereal*, doublereal*, void**);
	int (*psolve) (doublereal*, doublereal*, void**);
	void** vPData; //precond data
	void** vMData; //matvec data
} GMRESParams;

//The Params struct describes the problem setup, and administrivia
typedef struct
{
	//Administrative stuff
	double dT; //the total runtime
	int iM; //The number of timesteps
	double dTs; //time step size
	char sCase[20]; //I'm not sure if we will actually need this
	int iOrder; //The solver order
	FILE *fOut; //Output file
	//Physical stuff
	//something viscCont, I'm not sure what type this is...
	int bIncomp; //Toggles whether fluid is incompressible
	//something vInf, determines the far field velocity (in some functional way perhaps)
	double dLp;
} Params;

//The SingQuad struct contains the data needed to integrate a singularity
typedef struct
{
	double *dA;
	double *dWt;
	int iY; //rows of A and dWt
	int iX; //cols of A
} SingQuad;

//contains data to describe a velocity field
typedef struct
{
	char sFieldFun[20]; //ID for the velocity field function
	int NP; //number of parameters
	double *dMisc; //Stores misc parameters for the FieldFun
} VelField;

//Stores returns from GetBDF function
typedef struct
{
	double dBDF;
	double *dX; //This is for the X Coefficient
	double *dRHS; //This is for the RHS Coeff
} Coeffs;

//Stores returns from GetWeights function
typedef struct
{
	double *dV;
	double *dU;
	int iA;
} Weights;

//This is the return from the UpdateFirstM
//It's not really explained in the MATLAB code
typedef struct
{
	double *dX;
	double *dY;
	double *dSig;
	double *dAvgS;
} Store;

//stores the forces on an n point body
typedef struct
{
	int N;
	double *dFX;
	double *dFY;
} Force;

//This stores a vesicle and a problem case (velocity field, storage, forces, kernels, and SigM)
typedef struct
{
	int N; //number of vesicles
	Ves2D V; //the actual vesicle
	Store STORE; //records stuff from the last time step
	Force FS; //force from self interaction
	Force FO; //force from interaction with other vesicles
	Force FFar; //force from far field velocity
	Force FOut; //force that comes back from matvec
	Force FV; //intermediate force from interaction force
	Force FD; //double layer force
	double *dG; //G Kernel
	double *dD; //D Kernel
	double *dSigM; //I'm not really sure what this does
	double *dLambda;
	double *dL;
	double *dRHSVec; // this was pulled off of Coeffs
	//Parameters for GMRES
	GMRESParams GMParLHS;
	GMRESParams GMParTime;
	SingQuad SQ;//for quadrature around singularities
} VesCase;

//Convenient top level struct
typedef struct
{
	int NV;
	VesCase *VC;
} VesCases;

//handy type for passing around EvalVel
typedef void (*EvalVel)(int, double*, double*, double*, double*, int, ...);

//Maker functions
void Make_VesCase(int N, char *sDescriptor, VesCase *VC);
void Make_Force(int N, Force *F);
void Make_Store(int N, Store *STORE);

//These are used by Make_VesCase
double *GetLambda(int N);
double *GetL(int N);
void GetWeights(int iQ, int iKer, Weights* W);
void GetRegWeights(int iA, double *dX, double *dW);
void QuadLogSingData(int iM, int iQ, double dL, SingQuad *SQ);

//The rest
void EvalFFVelocity(int N, double *dX, double *dY, double *dU, double *dV, int iMode, ...);

//destroying my funny structs
void Dest_SingQuad(SingQuad *SQ);
void Dest_Coeffs(Coeffs *Co);
void Dest_Weights(Weights *W);
void Dest_Store(Store *STORE);
void Dest_Force(Force *F);
void Dest_VesCase(VesCase *VC);
void Dest_VesCases(VesCases *VCs);
void Dest_VelField(VelField *VF);
#endif
