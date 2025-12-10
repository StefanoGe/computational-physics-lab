#include "matrix.h"
#include "comp_physics.h"

void mat_init(Matrix *mat, int nrows, int ncols)
{
	mat_free(mat);
	mat->data = malloc(sizeof(double) * nrows * ncols);
	if(!mat->data) raiseErr("Error in memory allocation");
	for(int i=0; i<nrows; i++)
		mat->rows[i]=mat->data+i*ncols;
	mat->nrows = nrows;
	mat->ncols = ncols;
	mat->owns_data=true;
}

void mat_free( Matrix *mat )
{	
	if(mat->owns_data)
		free(mat->data);
	mat->data=nullptr;
	for(int i=0; i<mat->nrows; i++)
		mat->rows[i]=nullptr;
	mat->nrows=0;
	mat->ncols=0;
}

