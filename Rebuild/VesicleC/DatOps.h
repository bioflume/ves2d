/*******************************************************************************
	DatOps.h															6/21/09
	Alexander S Rattner

This has the functions to read and write dat files
*******************************************************************************/

#ifndef DATOPTS_H
#define DATOPTS_H

//Bring in Externals
#include <stdio.h>
#include <stdlib.h>

//Function prototypes
double **Make2D_double(int iY, int iX);
void Dest2D_double(double **D, int iY);
void ReadDat(char *sFile, double **dOut, int iLy, int iLx);

#endif
