#include "interpolation.h"
#include "linearsys.h"

int interp_polyn_vand_lup(const Array *x, const Array *y, Array *coeffs, 
	double *cond_number)
{
	if(x->size != y->size)
		raiseErr("x and y must be same size: x = %d y = %d", x->size, y->size);
		
	const int size = x->size;
	
	Matrix vandermont = mat_new(size,size);
	
	arr_init(coeffs,size);
	for(int i=0; i<coeffs->size; i++)
		ARRP(coeffs,i)=ARRP(y,i);
	
	// Set first col to 1 ( as it is the col of the constant function)
	for(int row = 0; row < size; row++ )
		MAT(vandermont,row,0) = 1.0;
	
	// Compute rest of Vandermonde matrix
	for( int col = 1; col < size; col++ )
		for( int row = 0; row < size; row++ )
			MAT(vandermont,row,col) = MAT(vandermont,row,col - 1) * ARRP(x,row);
	
	int *pivots =malloc(sizeof(int)*size);
	
	if(cond_number)
	{
		*cond_number = linst_cond_numb_inf(&vandermont);
	}
	
	if( linst_lup_factor(&vandermont, pivots, -1) != 0 )
		return SINGULAR;
	
	linst_lup_solve_inplace(&vandermont, pivots, coeffs);

	mat_free(&vandermont);
	free(pivots);
	return NONSINGULAR;
}
