#ifndef MATRIX_H
#define MATRIX_H

#define MAT(A,r,c) ((A).rows[(r)][(c)])
#define MATP(A,r,c) ((A)->rows[(r)][(c)])

typedef struct {
	double *data;
	double **rows;
	int nrows;
	int ncols;
	bool owns_data;
} Matrix;

void mat_init(Matrix *mat, int nrows, int ncols);

void mat_free( Matrix *mat );

void mat_build_rows(Matrix *mat);

Matrix mat_new(int nrows, int ncols);

void mat_transpose(Matrix *source, Matrix *dest);

//double *mat_getr_ptr(const Matrix *mat, int row);

#endif
