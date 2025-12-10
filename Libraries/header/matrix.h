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

double *mat_getr_ptr(const Matrix *mat, int row);
