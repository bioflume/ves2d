/*******************************************************************************
	Vesicle.h														6/16/09
	Alexander S Rattner

Defines the basic vesicle structure and geometry functions
*******************************************************************************/

#ifndef VESICLE_H
#define VESICLE_H

//Pull in externals
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <cblas.h>
#include <fftw3.h>
#include "Misc.h"
#include "Contour.h"

typedef struct
{
	int iID; //The current Vesicle ID
	char sDescriptor[MAX_DESC_LENGTH]; //I'm putting this in so that the vesicles can be named and distinguished 
	Contour C; //Contains all of the geometry data to define the boundary
	double dViscCont;
	double dKappa; //stiffness, not curvature
} Ves2D;

//Makers
void Make_Vesicle(char *sDescriptor, int N, Ves2D *V);

//Other functions
void PrintVesicle(Ves2D *V);

//Destroyers
void Dest_Vesicle(Ves2D *V);

#endif
