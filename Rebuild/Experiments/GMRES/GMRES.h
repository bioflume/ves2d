#ifndef GMRES_H
#define GMRES_H

#include "f2c.h"
#include <cblas.h>

int gmres(integer *n, doublereal *b, doublereal *x, integer *restrt, doublereal *work, integer *ldw, doublereal *h, integer *ldh, integer *iter, doublereal *resid, int (*matvec) (doublereal*, doublereal*, doublereal*, doublereal*, int, ...), int (*psolve) (doublereal*, doublereal*, int, ...), integer *info);

#endif