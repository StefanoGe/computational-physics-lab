#include "matrix.h"
#include "genutil.h"
#include <string.h>


void mat_init(Matrix *mat, int nrows, int ncols)
{
    if (mat->owns_data) {
        free(mat->data);
    }
	free(mat->rows);
	
	mat->data = malloc(sizeof(double) * nrows * ncols);
	if(!mat->data) raiseErr("Error in memory allocation");
	
	mat->nrows = nrows;
	mat->ncols = ncols;
	mat->owns_data=true;
	mat_build_rows(mat);
}

void mat_free( Matrix *mat )
{	
	free(mat->rows);
	if(mat->owns_data)
	{
		free(mat->data);
	}
	mat->nrows=0;
	mat->ncols=0;
	mat->data=nullptr;
	mat->rows=nullptr;
	mat->owns_data=false;
}

void mat_build_rows(Matrix *mat)
{
	mat->rows=malloc(sizeof(double*)*mat->nrows);
	if(!mat->rows) raiseErr("Error in memory allocation");
	for(int i=0; i<mat->nrows; i++)
		mat->rows[i]=mat->data+i*mat->ncols;
}

void mat_transpose(Matrix *source, Matrix *dest)
{
	const int new_nrows=source->ncols;
	const int new_ncols=source->nrows;
	double *new_data_ptr=malloc(sizeof(double)*new_ncols*new_nrows);
	if (!new_data_ptr) raiseErr("Error in memory allocation");
	
	for(int i=0; i<new_nrows; i++)
		for(int j=0; j<new_ncols; j++)
			new_data_ptr[i*new_ncols+j]=MATP(source,j,i);
	

	mat_free(dest);
	dest->data=new_data_ptr;
	dest->nrows=new_nrows;
	dest->ncols=new_ncols;
	dest->owns_data=true;
	mat_build_rows(dest);
}

Matrix mat_new(int nrows, int ncols)
{
	Matrix m={0};
	mat_init(&m, nrows, ncols);
	return m;
}

void mat_free_many(Matrix **mats, int n)
{
    for (int i = 0; i < n; i++)
        mat_free(mats[i]);
}

void mat_diag(Matrix *mat, double x)
{
	if(mat->ncols != mat->nrows)
		raiseErr("matrix is not nxn - nrows = %d ncols = %d", mat->nrows, mat->ncols);
	
	for(int row=0; row<mat->nrows; row++)
		for(int col=0; col<mat->ncols; col++)
			MATP(mat, row, col)= (row==col) ? x : 0;
}

Matrix mat_new_from_file(const char * filename, int nrows, int ncols)
{
	FILE * file = openFile( filename, "r" );
	
	Matrix mat = mat_new( nrows, ncols );
	for( int i = 0; i < nrows; i ++ )
		for( int j = 0; j < ncols; j++ )
			if (fscanf(file, "%lf", &MAT(mat, i, j)) != 1)
				raiseErr("Invalid or incomplete matrix file");
	
	fclose(file);
	return mat;
}

void mat_print_stdout( const Matrix *mat, char *format, bool newline)
{
	const int last_col=mat->ncols-1;
	const int last_row=mat->nrows-1;
	for(int row=0; row<mat->nrows; row++)
	{
		fputs( (row==0) ? "[ " : "| ", stdout );
		for(int col=0; col<mat->ncols; col++ )
		{
			printf(format, MATP(mat,row, col));
			
			if( col != last_col )
				fputs(", ", stdout);
		}
		fputs( (row != last_row) ? "\n" : " ]", stdout );
	}
	if(newline)
		putchar('\n');
}

bool mat_is_squared(const Matrix *m) {return m->ncols == m->nrows;}

void mat_swap_row(Matrix *m, int i, int j)
{
	if(i==j)
		return;
	double temp=0;
	for(int p=0; p<m->ncols; p++)
	{
		temp=MATP(m,i,p);
		MATP(m,i,p)=MATP(m,j,p);
		MATP(m,j,p)=temp;
	}
}

void mat_setv(const Matrix *A, double x)
{
	for(int row=0; row<A->nrows; row++)
		for(int col=0; col<A->ncols; col++)
			MATP(A, row, col)= x;
}

void mat_mult(const Matrix *A, const Matrix *B, Matrix *dest)
{	
	if( A->ncols!=B->nrows )
		raiseErr("A and B are not compatible for multiplication");
	
    bool is_dest_factor = (dest == A) || (dest == B);

    Matrix temp = {0};
    Matrix *out = dest;

    if (is_dest_factor)
    {
        mat_init(&temp, A->nrows, B->ncols);
        out = &temp;
    }
    else
    {
        mat_init(dest, A->nrows, A->ncols);
    }
	
	mat_setv(out,0);
	
	for(int i=0; i<A->nrows; i++)
		for(int j=0; j<B->ncols; j++)
			for(int p=0; p<A->ncols; p++)
				MATP(out,i,j)+=MATP(A,i,p)*MATP(B,p,j);
				
	if(is_dest_factor)
	{
		mat_cp(out,dest);
		mat_free(out);
	}
}

Matrix mat_mult_new(const Matrix *A, const Matrix *B)
{
	Matrix C=mat_new(A->nrows, B->ncols);
	mat_mult(A,B,&C);
	return C;	
}

void mat_vecmult(const Matrix *A, const Array *x, Array *dest)
{
	arr_init(dest, A->nrows);
	if(A->ncols!=x->size)
		raiseErr("A and B are not compatible for multiplication");

	arr_setv(dest, 0);
	for(int i=0; i<A->nrows; i++)
		for(int p=0; p<A->ncols; p++)
			ARRP(dest,i)+=MATP(A,i,p)*ARRP(x,p);
}

Array mat_vecmult_new(const Matrix *A, const Array *x)
{
	Array dest = arr_new(A->nrows);
	mat_vecmult(A,x,&dest);
	return dest;
}

void mat_cp(const Matrix *src, Matrix *dest)
{
	mat_init(dest,src->nrows, src->ncols);
	memcpy(dest->data, src->data, src->nrows*src->ncols*sizeof(double));
}

void mat_scale(Matrix *A, double alpha)
{
	for(int i=0; i<A->nrows*A->ncols; i++)
		A->data[i]*=alpha;
}

void mat_axpy(Matrix *y, double alpha, const Matrix *x)
{
	if(y->nrows!=x->nrows || y->ncols!=x->ncols)
		raiseErr("matrices should be of same dimension");
	for(int i=0; i<y->nrows*y->ncols; i++)
		y->data[i]+=alpha*x->data[i];
}

/*
 * Given a matrix A, constructs the matrix AtA (A transpose A, where first A
 * is transposed)
 * 
 * Args:
 * A: source
 * G: destination
 * 
 */

void mat_gram(const Matrix *A, Matrix *G)
{
	mat_init(G, A->ncols, A->ncols);
	
	for(int row=0; row< A->nrows; row++)
		for(int col=row; col<A->ncols; col++)
		{
			double sum=0;
			
			for(int p=0; p<A->nrows; p++)
				sum+=MATP(A,p,row)*MATP(A,p,col);
			
			MATP(G,row,col)=MATP(G,col,row)=sum;
		}
}

double mat_dot_col_col(const Matrix *m, int c1, int c2)
{
	double sum=0;
	for(int r=0; r<m->nrows; r++)
		sum+=MATP(m,r,c1)*MATP(m,r,c2);
	
	return sum;
}

double mat_dot_col_arr(const Matrix *m, int c, const Array *arr)
{
	double sum=0;
	for(int r=0; r<m->nrows; r++)
		sum+=MATP(m,r,c)*ARRP(arr,r);
	
	return sum;
}

void mat_aug_arr(const Matrix *A, const Array*b, Matrix *dest)
{
	if(A->nrows!=b->size)
		raiseErr("A and b must be of same height");
	mat_init(dest, A->nrows, A->ncols+1);
	for(int i=0; i<A->nrows; i++)
	{
		memcpy(*(dest->rows+i), *(A->rows+i), A->ncols*sizeof(double));
		MATP(dest,i,A->ncols)=ARRP(b,i);
	}
}

Matrix mat_aug_arr_new(const Matrix *A, const Array*b)
{
	Matrix dest={0};
	mat_aug_arr(A,b,&dest);
	return dest;
}

void mat_col_scale(Matrix *A, int c, double x)
{
	if(A->ncols<=c)
		raiseErr("c selects out of bounds column");
	for(int i=0; i<A->nrows; i++)
		MATP(A,i,c)=MATP(A,i,c)*x;
}

void mat_axpy_col_col(Matrix *y, int cy, double alpha, const Matrix *x, int cx)
{
	if(x->nrows!=y->nrows)
		raiseErr("mats should be of same height");
	for(int i=0; i<x->nrows; i++)
		MATP(y,i,cy)+=alpha*MATP(x,i,cx);
}

void mat_col_to_arr(const Matrix *mat, Array *arr, int col)
{
	if(col>=mat->ncols)
		raiseErr("Index out of bounds");
	arr_init(arr, mat->nrows);
	
	for(int i=0; i<mat->nrows; i++)
		ARRP(arr,i)=MATP(mat,i,col);
}


