/*******************************************************************************
	DatOps.c															6/21/09
	Alexander S Rattner

This has the functions to read and write dat files
*******************************************************************************/

#include "DatOps.h"

//Makes a 2D double array
double **Make2D_double(int iY, int iX)
{
	double **D;
	int i;
	
	D = malloc(sizeof(double *) * iY);
	for (i=0; i<iY; i++)
		D[i] = malloc(sizeof(double) * iX);
	
	return D;
}

//Destroys a 2D double array
void Dest2D_double(double **D, int iY)
{
	int i;
	
	for(i=0; i<iY; i++)
		free(D[i]);
	free(D);
}

//Reads a data file of iLx columns and iLy rows into dOut
void ReadDat(char *sFile, double **dOut, int iLy, int iLx)
{
	int i, j;
	FILE *fDat;
	
	fDat = fopen(sFile, "r");
	
	//Can't find file error handling
	if(fDat==NULL)
	{
		printf("Error: can't access dat file\n");
		return;
	}
	
	//read data into dOut
	for (i=0; i<iLy; i++)
		for (j=0; j<iLx; j++)
			fscanf(fDat, "%lf", &(dOut[i][j]));
	
	fclose(fDat);
}
