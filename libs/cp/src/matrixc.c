#include "matrixc.h"
#include "genutil.h"

/*
 * New matrices should always be 0-initialized. e.g.
 * Matrix mat={0};
 * 
 */

void matc_init(MatrixC *mat, int nrows, int ncols)
{
	matc_free(mat);
	
	mat->data = malloc(sizeof(double complex) * nrows * ncols);
	if(!mat->data) raiseErr("Error in memory allocation");
	
	mat->nrows = nrows;
	mat->ncols = ncols;
	mat->owns_data=true;
	matc_build_rows(mat);
}

void matc_free( MatrixC *mat )
{	
	if(!mat) return;
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

void matc_build_rows(MatrixC *mat)
{
	mat->rows=malloc(sizeof(double complex*)*mat->nrows);
	if(!mat->rows) raiseErr("Error in memory allocation");
	for(int i=0; i<mat->nrows; i++)
		mat->rows[i]=mat->data+i*mat->ncols;
}

MatrixC matc_new(int nrows, int ncols)
{
	MatrixC m={0};
	matc_init(&m, nrows, ncols);
	return m;
}

void matc_free_many(MatrixC **mats, int n)
{
    for (int i = 0; i < n; i++)
        matc_free(mats[i]);
}

void matc_set_value(MatrixC *mat, double complex value)
{
	for(int i=0; i<mat->nrows; i++)
		for(int j=0; j<mat->ncols; j++)
			MATP(mat,i,j)=value;
}

void matc_print(FILE *f, const MatrixC *mat)
{
	for(int row=0; row<mat->nrows; row++)
	{
		for(int col=0; col<mat->ncols; col++)
			fprintf(f, "%.2lf%+.2lfi  ", creal(MATP(mat,row,col)), cimag(MATP(mat,row,col)));
		fputs("\n", f);
	}
}

void matc_arrc_mult(double complex **mat, const double complex *arr, int nrows, 
				int ncols, ArrayC *dest)
{
	double complex *new_data=malloc(sizeof(double complex)*nrows);
	for(int i=0; i<nrows; i++)
	{
		new_data[i]=0;
		for(int j=0; j<ncols; j++)
			new_data[i]+=mat[i][j]*arr[j];
	}
	arrc_free(dest);
	dest->data=new_data;
	dest->owns_data=true;
	dest->size=nrows;
}







