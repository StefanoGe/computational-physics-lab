#include "matrix.h"
#include "comp_physics.h"

void mat_init(Matrix *mat, int nrows, int ncols)
{
	mat_free(mat);
	mat->data = malloc(sizeof(double) * nrows * ncols);
	if(!mat->data) raiseErr("Error in memory allocation");
	mat->rows = malloc(sizeof(double*)*nrows);
	for(int i=0; i<nrows; i++)
		mat->rows[i]=mat->data+i*ncols;
	mat->nrows = nrows;
	mat->ncols = ncols;
	mat->owns_data=true;
}

void mat_free( Matrix *mat )
{	
	if(mat->owns_data)
	{
		free(mat->data);
		free(mat->rows);
	}
	mat->nrows=0;
	mat->ncols=0;
}

