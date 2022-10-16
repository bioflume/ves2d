/*******************************************************************************
	TestTop.c														6/16/09
	Alexander S Rattner

Validates the performance of the Vesicle code.
*******************************************************************************/

//basic libraries
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <time.h>
#include <libxml/xmlmemory.h>
#include <libxml/parser.h>
#include <fftw3.h>

//custom includes
#include "Vesicle.h"
#include "Stokes.h"
#include "StokesStructs.h"
#include "Boundary.h"
#include "VesXML.h"

//Top level function
int main(int argc, char *argv[])
{
	Params Par; //solution parameters
	VesCases VCs; //top level problem cases
	Boundary B, *Bptr; //the boundary structure (but of course)
	VelField VF; //far field velocity (used in unbounded case)
	char *sFileIn = argv[1];
	char *sFileOut = argv[2];
	EvalVel EV;
	int iRes = -10;
	time_t tTime; //stores elapsed solver time
	
	//read input file
	iRes = ReadCase(sFileIn, &VCs, &B, &VF, &Par, &EV);
	
	//check for unbounded (then set Bptr to null)
	Bptr = (B.NC != -1) ? &B : (Boundary*)NULL;		
	
	//initialize output
	WriteOutput("Start", 3, Bptr, &VF, &VCs);
	
	printf("Stating solve...\n");
	tTime = time(NULL);
	
	//run study
	SOLVE(&VCs, &Par, Bptr, EV);
	
	printf("Solve done, elapsed time: %ds\n", (int)difftime(time(NULL), tTime));
	
	//make output file
	WriteOutput("End", 1, sFileOut);
	
	//clear up stuff
	Dest_VesCases(&VCs);
	if (B.NC != -1)
	{	Dest_Boundary(&B); }
	else
	{	Dest_VelField(&VF); }
	fftw_cleanup();
	xmlCleanupParser(); 
	
	return 0;
}
