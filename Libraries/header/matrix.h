#ifndef MATRIX_H
#define MATRIX_H

#ifndef MAT_ACCESS
#define MAT_ACCESS

#define MAT(A,r,c) ((A).rows[(r)][(c)])
#define MATP(A,r,c) ((A)->rows[(r)][(c)])

#endif

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

void mat_free_many(Matrix **mats, int n);

void mat_diag(Matrix *, double x);

Matrix mat_new_from_file(char * filename, int nrows, int ncols);

void mat_print_stdout( const Matrix *mat, char *format, bool newline);

#endif
