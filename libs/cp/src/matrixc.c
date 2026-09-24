#include "matrixc.h"
#include "genutil.h"
#include <string.h>


void matc_init(MatrixC *mat, int nrows, int ncols)
{
    if (mat->owns_data) {
        free(mat->data);
    }
    free(mat->rows);

    mat->data = malloc(sizeof(double complex) * nrows * ncols);
    if(!mat->data) raiseErr("Error in memory allocation");

    mat->nrows = nrows;
    mat->ncols = ncols;
    mat->owns_data = true;
    matc_build_rows(mat);
}

void matc_free(MatrixC *mat)
{
    free(mat->rows);
    if(mat->owns_data)
    {
        free(mat->data);
    }

    mat->nrows = 0;
    mat->ncols = 0;
    mat->data = NULL;
    mat->rows = NULL;
    mat->owns_data = false;
}

void matc_build_rows(MatrixC *mat)
{
    mat->rows = malloc(sizeof(double complex *) * mat->nrows);
    if(!mat->rows) raiseErr("Error in memory allocation");

    for(int i = 0; i < mat->nrows; i++)
        mat->rows[i] = mat->data + i * mat->ncols;
}

void matc_transpose(MatrixC *source, MatrixC *dest)
{
    const int new_nrows = source->ncols;
    const int new_ncols = source->nrows;

    double complex *new_data_ptr =
        malloc(sizeof(double complex) * new_ncols * new_nrows);

    if(!new_data_ptr)
        raiseErr("Error in memory allocation");

    for(int i = 0; i < new_nrows; i++)
        for(int j = 0; j < new_ncols; j++)
            new_data_ptr[i * new_ncols + j] = MATP(source, j, i);

    matc_free(dest);

    dest->data = new_data_ptr;
    dest->nrows = new_nrows;
    dest->ncols = new_ncols;
    dest->owns_data = true;

    matc_build_rows(dest);
}

MatrixC matc_new(int nrows, int ncols)
{
    MatrixC m = {0};
    matc_init(&m, nrows, ncols);
    return m;
}

void matc_free_many(MatrixC **mats, int n)
{
    for(int i = 0; i < n; i++)
        matc_free(mats[i]);
}

void matc_diag(MatrixC *mat, double complex x)
{
    if(mat->ncols != mat->nrows)
        raiseErr("matrix is not nxn - nrows = %d ncols = %d",
                 mat->nrows, mat->ncols);

    for(int row = 0; row < mat->nrows; row++)
        for(int col = 0; col < mat->ncols; col++)
            MATP(mat, row, col) = (row == col) ? x : 0;
}

MatrixC matc_new_from_file(const char *filename, int nrows, int ncols)
{
    FILE *file = openFile(filename, "r");

    MatrixC mat = matc_new(nrows, ncols);

    /*
     * Complex input format needs to be defined explicitly.
     * This currently assumes each value is stored as:
     *
     *     real imag
     *
     * for example:
     *
     *     1.0  2.0
     *     3.0 -1.0
     */
    for(int i = 0; i < nrows; i++)
        for(int j = 0; j < ncols; j++)
        {
            double re, im;

            if(fscanf(file, "%lf %lf", &re, &im) != 2)
                raiseErr("Invalid or incomplete matrix file");

            MAT(mat, i, j) = re + im * I;
        }

    fclose(file);
    return mat;
}

void matc_print_stdout(const MatrixC *mat, char *formatr,
	char *formati, bool newline)
{
    const int last_col = mat->ncols - 1;
    const int last_row = mat->nrows - 1;

    for(int row = 0; row < mat->nrows; row++)
    {
        fputs((row == 0) ? "[ " : "| ", stdout);

        for(int col = 0; col < mat->ncols; col++)
        {
            /*
             * printf cannot print a double complex directly.
             * The format is therefore applied separately to the
             * real and imaginary components.
             */
            printf(formatr, creal(MATP(mat, row, col)));
            printf(formati, cimag(MATP(mat, row, col)));

            if(col != last_col)
                fputs(", ", stdout);
        }

        fputs((row != last_row) ? "\n" : " ]", stdout);
    }

    if(newline)
        putchar('\n');
}

bool matc_is_squared(const MatrixC *m)
{
    return m->ncols == m->nrows;
}

void matc_swap_row(MatrixC *m, int i, int j)
{
    if(i == j)
        return;

    double complex temp = 0;

    for(int p = 0; p < m->ncols; p++)
    {
        temp = MATP(m, i, p);
        MATP(m, i, p) = MATP(m, j, p);
        MATP(m, j, p) = temp;
    }
}

void matc_setv(const MatrixC *A, double complex x)
{
    for(int row = 0; row < A->nrows; row++)
        for(int col = 0; col < A->ncols; col++)
            MATP(A, row, col) = x;
}

void matc_mult(const MatrixC *A, const MatrixC *B, MatrixC *dest)
{
    if(A->ncols != B->nrows)
        raiseErr("A and B are not compatible for multiplication");

    bool is_dest_factor = (dest == A) || (dest == B);

    MatrixC temp = {0};
    MatrixC *out = dest;

    if(is_dest_factor)
    {
        temp = matc_new(A->nrows, B->ncols);
        out = &temp;
    }
    else
    {
        matc_init(dest, A->nrows, B->ncols);
    }

    matc_setv(out, 0);

    for(int i = 0; i < A->nrows; i++)
        for(int j = 0; j < B->ncols; j++)
            for(int p = 0; p < A->ncols; p++)
                MATP(out, i, j) += MATP(A, i, p) * MATP(B, p, j);

    if(is_dest_factor)
    {
        matc_cp(out, dest);
        matc_free(out);
    }
}

MatrixC matc_mult_new(const MatrixC *A, const MatrixC *B)
{
    MatrixC C = matc_new(A->nrows, B->ncols);
    matc_mult(A, B, &C);
    return C;
}

void matc_vecmult(const MatrixC *A, const ArrayC *x, ArrayC *dest)
{
    if(A->ncols != x->size)
        raiseErr("A and x are not compatible for multiplication");

    arrc_init(dest, A->nrows);
    arrc_setv(dest, 0);

    for(int i = 0; i < A->nrows; i++)
        for(int p = 0; p < A->ncols; p++)
            ARRP(dest, i) += MATP(A, i, p) * ARRP(x, p);
}

ArrayC matc_vecmult_new(const MatrixC *A, const ArrayC *x)
{
    ArrayC dest = arrc_new(A->nrows);
    matc_vecmult(A, x, &dest);
    return dest;
}

void matc_cp(const MatrixC *src, MatrixC *dest)
{
    matc_init(dest, src->nrows, src->ncols);

    memcpy(dest->data, src->data,
           src->nrows * src->ncols * sizeof *dest->data);
}

void matc_scale(MatrixC *A, double complex alpha)
{
    for(int i = 0; i < A->nrows * A->ncols; i++)
        A->data[i] *= alpha;
}

void matc_axpy(MatrixC *y, double complex alpha, const MatrixC *x)
{
    if(y->nrows != x->nrows || y->ncols != x->ncols)
        raiseErr("matrices should be of same dimension");

    for(int i = 0; i < y->nrows * y->ncols; i++)
        y->data[i] += alpha * x->data[i];
}

/*
 * Given a matrix A, constructs the Hermitian Gram matrix A^H A.
 *
 * G(i,j) = sum_p conj(A(p,i)) * A(p,j)
 */

void matc_gram(const MatrixC *A, MatrixC *G)
{
    matc_init(G, A->ncols, A->ncols);

    for(int row = 0; row < A->ncols; row++)
        for(int col = row; col < A->ncols; col++)
        {
            double complex sum = 0;

            for(int p = 0; p < A->nrows; p++)
                sum += conj(MATP(A, p, row)) * MATP(A, p, col);

            MATP(G, row, col) = sum;
            MATP(G, col, row) = conj(sum);
        }
}

double complex matc_dot_col_col(const MatrixC *m, int c1, int c2)
{
    double complex sum = 0;

    for(int r = 0; r < m->nrows; r++)
        sum += conj(MATP(m, r, c1)) * MATP(m, r, c2);

    return sum;
}

double complex matc_dot_col_arr(const MatrixC *m, int c, const ArrayC *arr)
{
    double complex sum = 0;

    for(int r = 0; r < m->nrows; r++)
        sum += conj(MATP(m, r, c)) * ARRP(arr, r);

    return sum;
}

void matc_aug_arr(const MatrixC *A, const ArrayC *b, MatrixC *dest)
{
    if(A->nrows != b->size)
        raiseErr("A and b must be of same height");

    matc_init(dest, A->nrows, A->ncols + 1);

    for(int i = 0; i < A->nrows; i++)
    {
        memcpy(*(dest->rows + i),
               *(A->rows + i),
               A->ncols * sizeof *dest->data);

        MATP(dest, i, A->ncols) = ARRP(b, i);
    }
}

MatrixC matc_aug_arr_new(const MatrixC *A, const ArrayC *b)
{
    MatrixC dest = {0};
    matc_aug_arr(A, b, &dest);
    return dest;
}

void matc_col_scale(MatrixC *A, int c, double complex x)
{
    if(A->ncols <= c)
        raiseErr("c selects out of bounds column");

    for(int i = 0; i < A->nrows; i++)
        MATP(A, i, c) *= x;
}

void matc_axpy_col_col(MatrixC *y, int cy, double complex alpha,
                       const MatrixC *x, int cx)
{
    if(x->nrows != y->nrows)
        raiseErr("mats should be of same height");

    for(int i = 0; i < x->nrows; i++)
        MATP(y, i, cy) += alpha * MATP(x, i, cx);
}

void matc_col_to_arr(const MatrixC *mat, ArrayC *arr, int col)
{
    if(col >= mat->ncols)
        raiseErr("Index out of bounds");

    arrc_init(arr, mat->nrows);

    for(int i = 0; i < mat->nrows; i++)
        ARRP(arr, i) = MATP(mat, i, col);
}

void matc_carrc_mult(const MatrixC *m, const double complex *arr, ArrayC *dest)
{
	if(dest->size!=m->nrows||dest->data==NULL)
		arrc_init(dest, m->nrows);
	
	for(int i=0; i<m->nrows; i++)
	{
		ARRP(dest,i)=0;
		for(int j=0; j<m->ncols; j++)
			ARRP(dest,i)+=MATP(m,i,j)*arr[j];
	}
}

void matc_carrc_mult_banded(const MatrixC *m, const double complex *arr, 
	ArrayC *dest,int bandwidth)
{	
	if(!matc_is_squared(m))
		raiseErr("mat must be squared");
    
	if(dest->size!=m->nrows || dest->data == NULL)
		arrc_init(dest, m->nrows);
	
	int first_occupied_col=0;
	int last_occupied_col =0;
	
	for(int i=0; i<m->nrows; i++)
	{
		first_occupied_col=MAX(0,i-bandwidth);
		last_occupied_col=MIN(m->nrows-1,i+bandwidth);
		
		for(int j=first_occupied_col; j<=last_occupied_col; j++)
		{
			ARRP(dest,i)+=MATP(m,i,j)*arr[j];
		}
	}
}

#define BW 2
#define BW2 4
void matc_mult_banded2(const MatrixC *A, const MatrixC *B,
	MatrixC *dest)
{	
	if(!matc_is_squared(A)||!matc_is_squared(B))
		raiseErr("A and B must be squared");
		
	if(A->nrows!=B->nrows)
		raiseErr("A and B must have same dimension");
    
	if(dest->nrows!=A->nrows || dest->ncols != A->ncols || dest->data == NULL)
		matc_init(dest, A->nrows,B->ncols);

	matc_setv(dest,0);
	
	for(int i=0; i<A->nrows; i++)
	{
		
		for(int j=MAX(0,i-BW2); j<=MIN(B->ncols-1,i+BW2); j++)
		{
			const int pmin = MAX(0, MAX(i - BW, j - BW));
			const int pmax = MIN(A->nrows - 1, MIN(i + BW, j + BW));

			double complex sum = 0.0;

			for (int p = pmin; p <= pmax; ++p)
				sum += MATP(A, i, p) * MATP(B, p, j);

			MATP(dest, i, j) = sum;
		}
	}
}
#undef BW
#undef BW2

void matc_mult_banded(const MatrixC *A, const MatrixC *B,
	MatrixC *dest,int bandwidth)
{	
	if(!matc_is_squared(A)||!matc_is_squared(B))
		raiseErr("A and B must be squared");
		
	if(A->nrows!=B->nrows)
		raiseErr("A and B must have same dimension");
    
	if(dest->nrows!=A->nrows || dest->ncols != A->ncols || dest->data == NULL)
		matc_init(dest, A->nrows,A->ncols);
	
	matc_setv(dest,0);
	
	const int bandwidth2 = 2*bandwidth;
	
	for(int i=0; i<A->nrows; i++)
	{	
		for(int j=MAX(0,i-bandwidth2); j<=MAX(B->ncols,i-bandwidth2); j++)
		{
			for(int p=MAX(MAX(i-bandwidth2,j-bandwidth2),0); 
					p<MAX(A->nrows,MAX(i+bandwidth2,i+bandwidth2)); p++)
				MATP(dest,i,j)+=MATP(A,i,p)*MATP(B,p,j);
		}
	}
}

void arrc_cp_from_carr(const double complex *src, int size, ArrayC *dest)
{
	if( !dest->owns_data || dest->size!=size )
		arrc_init(dest, size);
	memcpy(dest->data, src, size*sizeof(double complex));
}

void matc_print_generic(const MatrixC *m,
	FILE *file,ComplexPrinter print_function,
	const char *format)
{
	for(int i=0; i<m->nrows; i++)
	{
		for(int j=0; j<m->ncols; j++)
			print_function(file,format,creal(MATP(m,i,j)),cimag(MATP(m,i,j)));
		fputc('\n',file);
	}
}

void matc_print_abs(const MatrixC *m, FILE *file,const char *format)
{
	for(int i=0; i<m->nrows; i++)
	{
		for(int j=0; j<m->ncols; j++)
			fprintf(file,format,cabs(MATP(m,i,j)));
		fputc('\n',file);
	}
}

Matrix matc_cabs(const MatrixC *m)
{
	Matrix new = mat_new(m->nrows,m->ncols);
	for(int i=0; i<m->nrows; i++)
		for(int j=0; j<m->ncols; j++)
			MAT(new,i,j)=cabs(MATP(m,i,j));
	
	return new;
}

Matrix matc_real_part(const MatrixC *m)
{
	Matrix new = mat_new(m->nrows,m->ncols);
	for(int i=0; i<m->nrows; i++)
		for(int j=0; j<m->ncols; j++)
			MAT(new,i,j)=creal(MATP(m,i,j));
	return new;
}

Matrix matc_imag_part(const MatrixC *m)
{
	Matrix new = mat_new(m->nrows,m->ncols);
	for(int i=0; i<m->nrows; i++)
		for(int j=0; j<m->ncols; j++)
			MAT(new,i,j)=cimag(MATP(m,i,j));
	return new;
}
