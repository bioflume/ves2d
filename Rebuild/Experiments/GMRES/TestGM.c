#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include "GMRES.h"
#include "f2c.h"

//Function prototypes
doublereal* GetB(int n);
double* GetA(int n);
void PrintVec(int n, doublereal* v);
int matvec(doublereal *alpha, doublereal *x, doublereal *beta, doublereal *y, int iMode, ...);
int psolve(doublereal *x, doublereal *b, int iMode, ...);

//int gmres(integer *n, doublereal *b, doublereal *x, integer *restrt, doublereal *work, integer *ldw, doublereal *h, 
//          integer *ldh, integer *iter, doublereal *resid, int (*matvec) (), int (*psolve) (), integer *info);

int main()
{
	//initializers
	size_t SZ;
	integer n;; //problem size
	double *A; //left hand side
	doublereal *b; //right hand side
	doublereal *x; //solution vec
	integer restrt; //restart value
	doublereal *work; //some workspace array
	integer ldw; //I don't know how big this should be
	doublereal *h; //hessenberg matrix
	integer ldh; //I don't know how big this should be
	integer iter; //number of iterations
	doublereal resid; //residual limit
	integer info; //stores some output
	
	//temp
	doublereal alpha;
	doublereal beta;
	alpha = 1;
	beta = 0;
	
	//fill in params
	n = 10;
	restrt = 7;
	ldw = 20;
	ldh = 10;
	iter = 200;
	resid = 1E-7;
	info = 0;
	
	//load in mat and vec
	b = GetB((int)n);
	A = GetA((int)n);
	
	//make room for x
	x = (doublereal*) calloc((size_t)n, sizeof(doublereal));
	
	//make room for work
	SZ = (size_t)(ldw * (restrt+4)) * sizeof(doublereal);
	work = (doublereal*) malloc(SZ);
	
	//make room for h
	SZ = (size_t)(ldh * (restrt+2)) * sizeof(doublereal);
	h = (doublereal*) malloc(SZ);
	
	//initialize matvec and psolve
	matvec((doublereal*)NULL, (doublereal*)NULL, (doublereal*)NULL, (doublereal*)NULL, 2, n, A);
	psolve((doublereal*)NULL, (doublereal*)NULL, 1, n);		
	
	//DO IT
	gmres(&n, b, x, &restrt, work, &ldw, h, &ldh, &iter, &resid, &matvec, &psolve, &info);
	
	PrintVec((int)n, x);	
	
	printf("%lf\n", resid);
	printf("%d\n", (int)iter);
	
	//free up space
	free(b);
	free(x);
	free(A);
	free(work);
	free(h);
	
	return 0;
}

//Function defs
doublereal* GetB(int n)
{
	doublereal *b;
	size_t SZ;
	int i;
	FILE *fDat;
	
	//make space for b
	SZ = (size_t)n * sizeof(doublereal);
	b = (doublereal*) malloc(SZ);
	
	//load file
	fDat = fopen("b.dat", "r");
	
	//fill in b
	for (i=0; i<n; i++)
		fscanf(fDat, "%lf", &( b[i] ) );
	
	fclose(fDat);
	return b;
}

double* GetA(int n)
{
	double *A;
	size_t SZ;
	int i;
	FILE *fDat;
	
	//make space for b
	SZ = (size_t)(n*n) * sizeof(double);
	A = (double*) malloc(SZ);
	
	//load file
	fDat = fopen("A.dat", "r");
	
	//fill in b
	for (i=0; i<(n*n); i++)
		fscanf(fDat, "%lf", &( A[i] ) );
	
	fclose(fDat);
	return A;
}

void PrintVec(int n, doublereal* v)
{
	int i;
	
	for (i=0; i<n; i++)
		printf("%1.4lf ", (double)v[i]);
	printf("\n");
}

int matvec(doublereal *alpha, doublereal *x, doublereal *beta, doublereal *y, int iMode, ...)
{
	va_list vaArgs;
	static int N;
	static double *dA;
	
	if (iMode == 2) //initialization mode
	{
		va_start (vaArgs, iMode);  
		N = (int) va_arg(vaArgs, integer);
		dA = va_arg(vaArgs, double*);
		va_end(vaArgs);
	}
	else
	{
		cblas_dgemv(CblasRowMajor, CblasNoTrans, N, N, (double)*alpha, dA, N, (double*)x, 1, (double)*beta, (double*)y, 1);
	}
	return 0;
}

int psolve(doublereal *x, doublereal *b, int iMode, ...)
{
	va_list vaArgs;
	static int N;
	size_t SZ;
	
	if (iMode == 1) //initialization mode
	{
		va_start (vaArgs, iMode);  
		N = (int) va_arg(vaArgs, integer);
		va_end(vaArgs);
	}
	else
	{
		SZ = (size_t)N * sizeof(doublereal);
		memcpy(x, b, SZ); //identity preconditioner
	}
	return 0;
}
