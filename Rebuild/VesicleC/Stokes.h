/*******************************************************************************
	Stokes.h															6/20/09
	Alexander S Rattner

Defines the basic elements needed to solve 2D Stokes flow
*******************************************************************************/

#ifndef STOKES_H
#define STOKES_H

//include externals here
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>
#include <cblas.h>
#include <clapack.h>
#include <string.h>
#include <fftw3.h>
#include "DatOps.h"
#include "Contour.h"
#include "Vesicle.h"
#include "Misc.h"
#include "StokesStructs.h"
#include "GMRES.h"
#include "f2c.h"
#include "VesXML.h"
#include "Boundary.h"

//Solver function definitions
void GetBDFcoeff(int iOrder, Coeffs *Co);
void PrecondL(int N, doublereal* b, doublereal* x, double* dLambda);
void PrecondM(int N, doublereal* b, doublereal* x, double* dL);
void UpdateSTORERHS(VesCase *VC, Coeffs *Co);
void UpdateFirstM(Ves2D *V, Params *Par, Store *STORE);
void ExactStokes(Force *F, double *dTX, double *dTY, int NT, double *dSX, double *dSY, int NS, Force *FV);
void InteractionForce(int NV, VesCase *VC);
void EvalDensity(Boundary *B, VesCases *VC);
void KernelS(Ves2D *V, SingQuad *SQ, double *dG);
void KernelD(Ves2D *V, double *dD);
void TimeStepper(Params *Par, Coeffs *Co, Force *FFar, double *dG, double *dD, Ves2D *V, double *dSig, 
                 Force *FOut, double *dRHSVec, GMRESParams *GMParLHS, GMRESParams *GMParTime);
double *SurfaceDiv(int N, double *dX, double *dTangX, double *dTangY, double *dS);
void LHSMatVec(double* dSigma, Ves2D *V, double *dG, double *dD, double *dOut);
void TimeMatVec(double *dX, Ves2D *V, double* dG, double *dD, Params *Par, Coeffs *Co, double *dVal, double *dSig, double *dFX, double *dFY, GMRESParams *GMParLHS);

//GMRES stuff
int LHSMatVecGM(doublereal *alpha, doublereal *x, doublereal *beta, doublereal *y, void** vMData);
int TimeMatVecGM(doublereal *alpha, doublereal *x, doublereal *beta, doublereal *y, void** vMData);
int PrecondLGM(doublereal *x, doublereal *b, void** vPData);
int PrecondMGM(doublereal *x, doublereal *b, void** vPData);
void RunGMRES(int iNin, double *dX, double *dB, GMRESParams *GMPar);

void AvgStress(Ves2D *V, Force *FTJ, double *dAvgS);

void SOLVE(VesCases *VCs, Params *Par, Boundary *B, EvalVel EV );
void SOLVE_VesStep(VesCase *VC, Params *Par, EvalVel EV, Coeffs *Co);

#endif
