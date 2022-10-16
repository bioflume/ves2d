/*******************************************************************************
	Vesicle.c														6/16/09
	Alexander S Rattner

Defines the basic vesicle structure and geometry functions
*******************************************************************************/

#include "Vesicle.h"

//Stick Globals here, is there a good way to avoid using this global??
int iCurID = 1;

//MAKER FUNCTIONS GO HERE

//Simple Vesicle Constructor
void Make_Vesicle(char *sDescriptor, int N, Ves2D *V)
{
	//Local Variable Defines
	int iNameLength;
	
	//Give it an ID
	V->iID = iCurID;
	iCurID++;
	
	//Give it a name, first check for valid length
	iNameLength = strlen(sDescriptor);
	if ( iNameLength <= MAX_DESC_LENGTH )   //acceptable name
	{
		strcpy(V->sDescriptor, sDescriptor);
	}
	else   //truncate the name
	{
		strncpy(V->sDescriptor, sDescriptor, MAX_DESC_LENGTH);
	}
	
	//Make the contour struct
	Make_Contour(N, &(V->C));
}

//MISC FUNCTIONS GO HERE

//Prints data about the vesicle
void PrintVesicle(Ves2D *V)
{
	int i;
	
	printf("\n============================================================\n");
	//Print out metadata
	printf("VesicleID: %d\t\t VesicleName: %s\n", V->iID, V->sDescriptor);
	printf("Number of Points: %d", V->C.N);
	printf("\t\tCenter Point: (%1.4lf, %1.4lf)\n", V->C.dCX, V->C.dCY);
	//Print out boundary data
	for (i=0; i<V->C.N; i++)
		printf("(%1.4lf, %1.4lf) ", V->C.dXPts[i], V->C.dYPts[i]);
	printf("\n============================================================\n");	
}

//DESTROYER FUNCTIONS GO HERE

//Clears Vesicle from Memory
void Dest_Vesicle(Ves2D *V)
{
	Dest_Contour( &(V->C) );
	// free(V);  This line won't be useful unless we malloced the Vesicle
}
