#include "matrix.h"
#include "genutil.h"

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

Matrix mat_new_from_file(char * filename, int nrows, int ncols)
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

