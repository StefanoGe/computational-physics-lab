#include "matrix.h"
#include "comp_physics.h"

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

