#include "linearsys.h"
#include "genutil.h"


void linst_forwsubst( double **ltm, int n, const double *b, Array *sol)
{	
	if( sol->size != n )
		arr_init(sol, n);
	
	double partial_sum;
	
	for( int row=0; row < n; row++ )
	{
		if (ltm[row][row] == 0.0)
			raiseErr("Singular matrix in forward substitution. row %d = %lf",
						row, ltm[row][row]);
		
		partial_sum = 0;
		for( int col = 0; col < row; col++ )
			partial_sum += ARRP(sol,col) * ltm[row][col];
		ARRP(sol,row) = (b[row] - partial_sum)/ltm[row][row];
	}
}

void linst_backsubst( double **utm, int n, const double *b, Array *sol)
{
	if( sol->size != n )
		arr_init(sol, n);
	
	double partial_sum;
	
	for( int row = n - 1; row >= 0; row-- )
	{
		if (utm[row][row] == 0.0)
			raiseErr("Singular matrix in forward substitution. row %d = %lf",
						row, utm[row][row]);
		
		partial_sum = 0;
		for( int col = row + 1; col < n; col++ )
			partial_sum += ARRP(sol,col) * utm[row][col];
		ARRP(sol,row) = (b[row] - partial_sum)/utm[row][row];
	}
}
