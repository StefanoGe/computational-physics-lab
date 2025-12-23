#ifndef MATRIXC_H
#define MATRIXC_H

#include <complex.h>
#include <stdio.h>
#include "arrayc.h"
#include <stdbool.h>

#ifndef MAT_ACCESS
#define MAT_ACCESS

	#define MAT(A,r,c) ((A).rows[(r)][(c)])
	#define MATP(A,r,c) ((A)->rows[(r)][(c)])

#endif


typedef struct {
	double complex *data;
	double complex **rows;
	int nrows;
	int ncols;
	bool owns_data;
} MatrixC;

void matc_init(MatrixC *mat, int nrows, int ncols);

void matc_free( MatrixC *mat );

void matc_build_rows(MatrixC *mat);

MatrixC matc_new(int nrows, int ncols);

void matc_free_many(MatrixC **mats, int n);

void matc_set_value(MatrixC *mat, double complex value);

void matc_print(FILE *f, const MatrixC *mat);

void matc_arrc_mult(double complex **mat, const double complex *arr, int nrows, 
				int ncols, ArrayC *dest);

#endif
