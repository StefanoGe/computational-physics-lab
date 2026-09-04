#ifndef MATRIXC_H
#define MATRIXC_H

#include <complex.h>
#include "arrayc.h"

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

void matc_free(MatrixC *mat);

void matc_build_rows(MatrixC *mat);

MatrixC matc_new(int nrows, int ncols);

void matc_transpose(MatrixC *source, MatrixC *dest);

void matc_free_many(MatrixC **mats, int n);

void matc_diag(MatrixC *, double complex x);

MatrixC matc_new_from_file(const char *filename, int nrows, int ncols);

void matc_print_stdout(const MatrixC *mat, char *formatr,
	char *formati, bool newline);

bool matc_is_squared(const MatrixC *m);

void matc_swap_row(MatrixC *m, int i, int j);

void matc_mult(const MatrixC *A, const MatrixC *B, MatrixC *dest);

MatrixC matc_mult_new(const MatrixC *A, const MatrixC *B);

void matc_vecmult(const MatrixC *A, const ArrayC *x, ArrayC *dest);

ArrayC matc_vecmult_new(const MatrixC *A, const ArrayC *x);

void matc_cp(const MatrixC *src, MatrixC *dest);

void matc_axpy(MatrixC *y, double complex alpha, const MatrixC *x);

void matc_scale(MatrixC *A, double complex alpha);

void matc_gram(const MatrixC *A, MatrixC *G);

void matc_setv(const MatrixC *A, double complex x);

double complex matc_dot_col_col(const MatrixC *m, int c1, int c2);

double complex matc_dot_col_arr(const MatrixC *m, int c, const ArrayC *arr);

void matc_aug_arr(const MatrixC *A, const ArrayC *b, MatrixC *dest);

MatrixC matc_aug_arr_new(const MatrixC *A, const ArrayC *b);

void matc_col_scale(MatrixC *A, int c, double complex x);

void matc_axpy_col_col(MatrixC *y, int cy, double complex alpha,
                       const MatrixC *x, int cx);

void matc_col_to_arr(const MatrixC *mat, ArrayC *arr, int col);

#endif
