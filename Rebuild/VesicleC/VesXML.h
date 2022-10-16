/*******************************************************************************
	VesXML.h															8/5/09
	Alexander S Rattner

Reads input files and writes output files
*******************************************************************************/

#ifndef VESXML_H
#define VESXML_H

//bring in externals
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <libxml/xmlmemory.h>
#include <libxml/parser.h>
#include <stdarg.h>
#include "Contour.h"
#include "Vesicle.h"
#include "Boundary.h"
#include "StokesStructs.h"
#include "Misc.h"
#include "Stokes.h"

//function definitions
int ReadCase(char* sFileIn, VesCases *VC, Boundary *B, VelField *VF, Params *Par, EvalVel *EV );
int WriteOutput(char* sMode, int iMode, ...);

//individual struct readers
void ReadVesCases(xmlNodePtr xVes, VesCases *VCs);
void ReadFarField(xmlNodePtr xFF, Boundary *B, VelField *VF, EvalVel *EV );
void ReadParams(xmlNodePtr xPar, Params *Par);

//individual struct writers
void WriteFarField(xmlNodePtr xFF, Boundary *B, VelField *VF);
void WriteVesCases(xmlNodePtr xVCs, VesCases *VCs);

//my handy xlib cheater functions
xmlNodePtr AddChild(xmlNodePtr xPar, char *sName);
void AddChildData(xmlNodePtr xPar, char *sName, char* sData);
void AddChildInt(xmlNodePtr xPar, char *sName, int iData);
void AddChildDouble(xmlNodePtr xPar, char *sName, double dData);
void AddChildDoubleArray(xmlNodePtr xPar, char *sName, double *dIn, int N);
xmlNodePtr GetChild(xmlNodePtr xPar, char* sName);
xmlChar *GetChildData(xmlNodePtr xPar, char* sName);
int GetChildInt(xmlNodePtr xPar, char* sName);
double GetChildDouble(xmlNodePtr xPar, char* sName);
void GetChildDoubleArray(xmlNodePtr xPar, char *sName, double *dOut, int N);
xmlNodePtr GetNeighbor(xmlNodePtr xOne, char* sName);

#endif
