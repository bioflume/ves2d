/*******************************************************************************
	VesXML.c															8/5/09
	Alexander S Rattner

Reads input files and writes output files
*******************************************************************************/

#include "VesXML.h"

//reads in test case xml file
int ReadCase(char* sFileIn, VesCases *VCs, Boundary *B, VelField *VF, Params *Par, EvalVel *EV )
{
	//XML Stuff
	xmlDocPtr xDoc;
	xmlNodePtr xRoot, xCur;
	
	printf("Reading input file: %s\n", sFileIn);
	
	//read in file and parse it
	xDoc = xmlParseFile(sFileIn);
	if (xDoc == NULL) //oops
	{
		printf("I can't read that input file...\n");
		return -1;
	}
	
	//start at root element (Study)
	xRoot = xmlDocGetRootElement(xDoc);
	if (xRoot == NULL)
	{
		printf("Empty doc\n");
		xmlFreeDoc(xDoc);
		return -1;
	}

	//read in ves cases
	xCur = GetChild(xRoot, "Vesicles");
	ReadVesCases(xCur, VCs);
	
	//Read in farfield stuff
	xCur = GetChild(xRoot, "FarField");
	ReadFarField(xCur, B, VF, EV);
	
	//Read in Params
	xCur = GetChild(xRoot, "Params");
	ReadParams(xCur, Par);
	
	//clean up, and return OK
	xmlFreeDoc(xDoc);
	return 0;
}

//Reads VesCases into VC
void ReadVesCases(xmlNodePtr xVes, VesCases *VCs)
{
	int NV;
	int N;
	int i;
	size_t SZ;
	xmlNodePtr xCGM, xCVes;
	
	//allocate vesicle memory
	NV = GetChildInt(xVes, "NV");
	VCs->NV = NV;
	SZ = (size_t)NV * sizeof(VesCase);
	VCs->VC = (VesCase*) malloc(SZ);
	
	//loop through vesicles
	xCVes = GetChild(xVes, "Vesicle");
	for (i=0; i<NV; i++)
	{
		//get number of points, and make VesCase
		N = GetChildInt(xCVes, "N");
		Make_VesCase(N, "", &(VCs->VC[i]));
		
		//read in Kappa
		VCs->VC[i].V.dKappa = GetChildDouble(xCVes, "Kappa");
		VCs->VC[i].V.dViscCont = 1.0;
		
		//read in X,Y data
		GetChildDoubleArray(xCVes, "X", VCs->VC[i].V.C.dXPts, N);
		GetChildDoubleArray(xCVes, "Y", VCs->VC[i].V.C.dYPts, N);
		
		//Get LHS params
		xCGM = GetChild(xCVes, "GMParLHS");
		VCs->VC[i].GMParLHS.iRestart = GetChildInt(xCGM, "Restart");
		VCs->VC[i].GMParLHS.iLdw = GetChildInt(xCGM, "LDW");
		VCs->VC[i].GMParLHS.iLdh = GetChildInt(xCGM, "LDH");
		VCs->VC[i].GMParLHS.iIterMax = GetChildInt(xCGM, "MaxIt");
		VCs->VC[i].GMParLHS.dResMax = GetChildDouble(xCGM, "MaxRes");
		VCs->VC[i].GMParLHS.matvec = &LHSMatVecGM;
		VCs->VC[i].GMParLHS.psolve = &PrecondLGM;
		
		//Get Time params
		xCGM = GetChild(xCVes, "GMParTime");
		VCs->VC[i].GMParTime.iRestart = GetChildInt(xCGM, "Restart");
		VCs->VC[i].GMParTime.iLdw = GetChildInt(xCGM, "LDW");
		VCs->VC[i].GMParTime.iLdh = GetChildInt(xCGM, "LDH");
		VCs->VC[i].GMParTime.iIterMax = GetChildInt(xCGM, "MaxIt");
		VCs->VC[i].GMParTime.dResMax = GetChildDouble(xCGM, "MaxRes");
		VCs->VC[i].GMParTime.matvec = &TimeMatVecGM;
		VCs->VC[i].GMParTime.psolve = &PrecondMGM;
		
		//Get next Vesicle
		if ( i < (NV-1) )
			xCVes = GetNeighbor(xCVes, "Vesicle");
	}
}

//Reads in the Far Field stuff
void ReadFarField(xmlNodePtr xFF, Boundary *B, VelField *VF, EvalVel *EV )
{
	xmlChar *xType, *xField;
	xmlNodePtr xCon, xPar;
	int NC, N, i;
	size_t SZ;
	
	//get type
	xType = GetChildData(xFF, "Type");
	
	if ( !xmlStrcmp(xType, (xmlChar*)"Bounded") ) //Bounded flow
	{
		//get Bounded top level node
		xFF = GetChild(xFF, "Boundary");
		
		//make the Contours
		NC = GetChildInt(xFF, "NC");
		SZ = (size_t)NC * sizeof(Contour);
		B->C = (Contour*) malloc(SZ);
		
		
		//fill in the contour data
		xCon = GetChild(xFF, "Contour");
		for (i=0; i<NC; i++)
		{
			//Initialize the Contour
			N = GetChildInt(xCon, "N");
			Make_Contour (N, &(B->C[i]));
			
			//read in the boundary data
			GetChildDoubleArray(xCon, "X", B->C[i].dXPts, N);
			GetChildDoubleArray(xCon, "Y", B->C[i].dYPts, N);
			GetChildDoubleArray(xCon, "U", B->C[i].dU, N);
			GetChildDoubleArray(xCon, "V", B->C[i].dV, N);
		
			//get next Contour
			if (i < (NC-1))
				xCon = GetNeighbor(xCon, "Contour");
		}
		
		//make the boundary
		Make_Boundary(1, NC, B->C, B);
		
		//Set the boundary function
		EvalBDVelocity((int)NULL, (double*)NULL, (double*)NULL, (double*)NULL, (double*)NULL, 1, B);
		*EV = &EvalBDVelocity;
	}
	else if ( !xmlStrcmp(xType, (xmlChar*)"Unbounded") ) //Unbounded flow
	{
		B->NC = -1; //signal Back for unbounded case
		xField = GetChildData(xFF, "FieldFun");
		
		//get flow type and data
		strcpy(VF->sFieldFun, (char*)xField);
		if ( !strcmp(VF->sFieldFun, "Shear") ) //shear flow, single parameter
		{	N = 1; }
		
		//initialize parameters array
		VF->NP = N;
		SZ = (size_t)N * sizeof(double);
		VF->dMisc = (double*) malloc(SZ);
		
		//read in parameters
		xPar = GetChild(xFF, "Param");
		for (i=0; i<N; i++)
		{
			sscanf( (char*)xmlNodeGetContent(xPar), "%lf", &(VF->dMisc[i]) );
			
			//get next parameter
			if (i < (N-1) )
				xPar = GetNeighbor(xPar, "Param");
		}
		
		//Initialize and set velocity function
		EvalFFVelocity((int)NULL, (double*)NULL, (double*)NULL, (double*)NULL, (double*)NULL, 1, VF);
		*EV = &EvalFFVelocity;
	
		//local clean up
		free(xField);
	}
	
	//clean up
	xmlFree(xType);
}

//Reads in params
void ReadParams(xmlNodePtr xPar, Params *Par)
{
	Par->dT = GetChildDouble(xPar, "Time");
	Par->iM = GetChildInt(xPar, "Steps");
	Par->dTs = Par->dT / (double)Par->iM;
	Par->iOrder = GetChildInt(xPar, "Order");
	Par->bIncomp = GetChildInt(xPar, "Incompressible");
	strcpy(Par->sCase, "schurComp");
	Par->dLp = 2*M_PI;
}

//writes output to file
int WriteOutput(char* sMode, int iMode, ...)
{
	va_list vaArgs;
	char *sFileOut;
	//stuff for making the file
	static xmlDocPtr xDoc;
	static xmlNodePtr xRoot;
	xmlNodePtr xFF, xVCs;
	double dT;
	
	//Data structs
	Boundary *B;
	VelField *VF;
	VesCases *VCs;
	
	//static xmlNodePtr xCur;
	
	if ( !strcmp(sMode, "Start") ) //Initialization mode
	{
		//initialize document
		xDoc = xmlNewDoc((const xmlChar*)"1.0");
		xRoot = xmlNewDocNode(xDoc, NULL, (const xmlChar*)"Results", NULL);
		xmlDocSetRootElement(xDoc, xRoot);
		
		//read in data
		va_start (vaArgs, iMode);  
		B = va_arg(vaArgs, Boundary*);
		VF = va_arg(vaArgs, VelField*);
		VCs = va_arg(vaArgs, VesCases*);
		va_end(vaArgs);
		
		//add Boundary data
		xFF = AddChild(xRoot, "FarField");
		WriteFarField(xFF, B, VF);
		
		//add VesCases data
		xVCs = AddChild(xRoot, "TimeStep");
		dT = 0.0; //initial time = 0
		AddChildDouble(xVCs, "T", dT); 
		WriteVesCases(xVCs, VCs);
	}
	else if ( !strcmp(sMode, "Step") ) //Add current time step 
	{
		//read in data
		va_start (vaArgs, iMode);
		VCs = va_arg(vaArgs, VesCases*);
		dT = va_arg(vaArgs, double);
		va_end(vaArgs);
		xVCs = AddChild(xRoot, "TimeStep");
		AddChildDouble(xVCs, "T", dT); 
		WriteVesCases(xVCs, VCs);
	}
	else if ( !strcmp(sMode, "End") ) //write file and close
	{
		//get output file name
		va_start (vaArgs, iMode);
		sFileOut = va_arg(vaArgs, char*);
		va_end(vaArgs);
		
		printf("Writing output file: %s\n", sFileOut);
		
		//write out
		xmlSaveFormatFile(sFileOut, xDoc, 1);
		//clean up
		xmlFreeDoc(xDoc);
	}
	
	return 0;
}

//writes far field data to XML struct
void WriteFarField(xmlNodePtr xFF, Boundary *B, VelField *VF)
{
	int i;
	int NC, N;
	xmlNodePtr xBnd, xCon;
	
	if (B != NULL) //Bounded flow
	{
		AddChildData(xFF, "Type", "Bounded");
		xBnd = AddChild(xFF, "Boundary");
		
		//basic boundary data
		NC = B->NC;
		AddChildInt(xBnd, "NC", NC);
		
		
		//start adding contours
		xCon = AddChild(xBnd, "Contour");
		for (i=0; i<NC; i++)
		{
			N = B->C[i].N;
			AddChildInt(xCon, "N", B->C[i].N);
			
			//add contour vectors
			AddChildDoubleArray(xCon, "X", B->C[i].dXPts, N);
			AddChildDoubleArray(xCon, "Y", B->C[i].dYPts, N);
			AddChildDoubleArray(xCon, "U", B->C[i].dU, N);
			AddChildDoubleArray(xCon, "V", B->C[i].dV, N);
			
			//start next Contour
			if ( i < (NC-1) )
				xCon = AddChild(xBnd, "Contour");
		}
	}
	else if (VF != NULL) //unbounded flow, print out the filed properties
	{
		N = VF->NP;
		AddChildData(xFF, "FieldFun", VF->sFieldFun);
		for (i=0; i<N; i++)
			AddChildDouble(xFF, "Param", VF->dMisc[i]);
	}
}

void WriteVesCases(xmlNodePtr xVCs, VesCases *VCs)
{
	int i = 0;
	int N, NV;
	xmlNodePtr xV;
	Contour *C;
	
	//print out vesicle data
	NV = VCs->NV;
	AddChildInt(xVCs, "NV", NV);
	xV = AddChild(xVCs, "Vesicle");
	for (i=0; i<NV; i++)
	{
		N = VCs->VC[i].V.C.N;
		AddChildInt(xV, "N", N);
		//print out position and velocity
		C = &(VCs->VC[i].V.C);
		AddChildDoubleArray(xV, "X", C->dXPts, N);
		AddChildDoubleArray(xV, "Y", C->dYPts, N);
		AddChildDoubleArray(xV, "U", C->dU, N);
		AddChildDoubleArray(xV, "V", C->dV, N);
		
		//start next Vesicle
		if ( i < (NV-1) )
			xV = AddChild(xVCs, "Vesicle");
	}
}

//adds a child element under the parent
xmlNodePtr AddChild(xmlNodePtr xPar, char *sName)
{
	xmlNodePtr xOut;
	xOut = xmlNewChild(xPar, NULL, (xmlChar*)sName, NULL);
	return xOut;
}

//adds a child element with a string value
void AddChildData(xmlNodePtr xPar, char *sName, char* sData)
{
	xmlNewChild(xPar, NULL, (xmlChar*)sName, (xmlChar*)sData );
}

//adds a child element with int value
void AddChildInt(xmlNodePtr xPar, char *sName, int iData)
{
	xmlChar xData[32]; //I hope that ints can't print longer than this
	sprintf( (char*)xData, "%d", iData);
	xmlNewChild(xPar, NULL, (xmlChar*)sName, xData);
}

//adds a child element with int value
void AddChildDouble(xmlNodePtr xPar, char *sName, double dData)
{
	xmlChar xData[32]; //I hope that doubles can't print longer than this
	sprintf( (char*)xData, "%lf", dData);
	xmlNewChild(xPar, NULL, (xmlChar*)sName, xData);
}

//adds a double array
void AddChildDoubleArray(xmlNodePtr xPar, char *sName, double *dIn, int N)
{
	int i, l, maxl, ind;
	int iDLength = 20; //I'm guessing this a reasonable length for printing doubles
	int iTLength; //the total length
	size_t SZ;
	xmlChar *xData;
	
	//initialize to first guess length
	iTLength = N*iDLength;
	SZ = (size_t)iTLength * sizeof(xmlChar);
	xData = (xmlChar*) xmlMalloc(SZ);
	
	//start sprintfing doubles
	ind=0; l=0; maxl = 0;
	for (i=0; i<N; i++) 
	{
		//write next double
		l = sprintf( (char*)&(xData[ind]), "%lf ", dIn[i] );
		ind += l;
		
		//check to see if xData is to small
		maxl = (maxl > l) ? maxl : l; //use the biggest length as a test size
		if ( maxl*(N-i) > (iTLength-ind) ) //we probably need more space
		{
			iTLength = maxl*N;
			SZ = (size_t)iTLength * sizeof(xmlChar);
			xData = realloc(xData, SZ);
		}
	}
	xmlNewChild(xPar, NULL, (xmlChar*)sName, xData);
	
	//clean up
	xmlFree(xData);
}

//returns child element under the parent element
xmlNodePtr GetChild(xmlNodePtr xPar, char* sName)
{
	xmlNodePtr xCur = xPar->children;
	
	while ((xCur != NULL) && xmlStrcmp(xCur->name, (xmlChar*) sName))
		xCur = xCur->next;
	
	return xCur;
}

//gets data from child
xmlChar *GetChildData(xmlNodePtr xPar, char* sName)
{
	xmlNodePtr xCur = GetChild(xPar, sName);
	return xmlNodeGetContent(xCur);
}

//gets int data from child
int GetChildInt(xmlNodePtr xPar, char* sName)
{
	xmlNodePtr xCur = GetChild(xPar, sName);
	xmlChar *xData = xmlNodeGetContent(xCur);
	int i;
	sscanf((char*)xData, "%d", &i);
	xmlFree(xData);
	return(i);
}

//gets int data from child
double GetChildDouble(xmlNodePtr xPar, char* sName)
{
	xmlNodePtr xCur = GetChild(xPar, sName);
	xmlChar *xData = xmlNodeGetContent(xCur);
	double d;
	sscanf((char*)xData, "%lf", &d);
	xmlFree(xData);
	return(d);
}

//reads child data into an array
void GetChildDoubleArray(xmlNodePtr xPar, char *sName, double *dOut, int N)
{
	xmlChar *xData = GetChildData(xPar, sName);
	ParseArray( (char*)xData, dOut, N);
	xmlFree(xData);
}

//gets proceding neighbor named sName
xmlNodePtr GetNeighbor(xmlNodePtr xOne, char* sName)
{
	xmlNodePtr xCur = xOne->next;
	
	while ((xCur != NULL) && xmlStrcmp(xCur->name, (xmlChar*) sName))
		xCur = xCur->next;
	
	return xCur;
}
