#ifndef LINEARSYS
#define LINEARSYS

#include "matrix.h"
#include "array.h"

void linst_forwsubst( double **lm, int n, const double *b, Array *sol);

void linst_backsubst( double **utm, int n, const double *b, Array *sol);

#endif
