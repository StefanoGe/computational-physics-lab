#define MAT(M,i,j) \
    (( (M)->data )[ (size_t)(i) * (M)->ncols + (j) ])

typedef struct {
	double *data;
	int nrows;
	int ncols;
} Matrix;

void mat_init(Matrix *mat, int nrows, int ncols);

void mat_free( Matrix *mat );

double *mat_getr_ptr(const Matrix *mat, int row);
