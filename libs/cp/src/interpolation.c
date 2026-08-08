#include "interpolation.h"
#include "linearsys.h"
#include <math.h>

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
	
	int *pivots;
	SAFE_ALLOC(pivots,size);
	
	if(cond_number)
	{
		*cond_number = linst_cond_numb_inf(&vandermont);
	}
	
	if( linst_lup_factor(&vandermont, pivots, -1) != 0 )
		return SINGULAR;
	
	linst_lup_solve_inplace(&vandermont, pivots, coeffs);

	mat_free(&vandermont);
	SAFE_FREE(pivots);
	return NONSINGULAR;
}

static inline double* interp_build_weights( const double *points, int npoints)
{
	double *weights;
	SAFE_ALLOC(weights,npoints);
	double prod = 1;

	for( int w = 0; w < npoints; w++ )
	{
		prod = 1;
		for( int j = 0; j < npoints; j++ )
			if( w != j )
				prod *= ( points[w] - points[j] );
				
		weights[w] = 1 / prod;
	}
	return weights;
}


void interp_barf_init_points( BarFit *barf , const Array *points)
{
	barf->npoints = points->size;
	SAFE_ALLOC(barf->points,barf->npoints);
	for(int i=0; i<barf->npoints; i++)
		barf->points[i]=ARRP(points,i);
	barf->weights = interp_build_weights( points->data, barf->npoints );
}

BarFit interp_barf_new_points(const Array *points)
{
	BarFit barf={0};
	interp_barf_init_points(&barf, points);
	return barf;
}

static inline double* interp_build_eq_weights(int n)
{
	double *weights;
	SAFE_ALLOC(weights,n);
	weights[ 0 ] = 1;
	for( int i = 1; i < n; i++ )
		weights[i] = weights[ i -1 ] * ( i - n ) / i;
	return weights;
}

void interp_barf_init_eq_points( BarFit *barf, double x1, double x2, int npoints )
{
	barf->npoints = npoints;
	SAFE_ALLOC(barf->points,npoints);
	const double step = (x2 - x1)/(npoints-1);
	for(int i=0; i<npoints; i++)
		barf->points[i]=x1 + step*i;
	barf->weights = interp_build_eq_weights( npoints );
}

BarFit interp_barf_new_eq_points(double x1, double x2, int npoints)
{
	BarFit barf={0};
	interp_barf_init_eq_points(&barf, x1, x2, npoints);
	return barf;
}

void interp_barf_add_func(BarFit *barf, const ParamFunc *f)
{
	if(!f)
		return;
	barf->func=f->func;
	barf->params=f->params;	
}

void interp_barf_add_fvalues(BarFit *barf, const ParamFunc *opt_f)
{	
	if(!barf->points||barf->npoints<=0)
		raiseErr("points are missing. cannot add fvalues");
		
	if(!barf->weights)
		raiseErr("weights are missing. cannot add fvalues");
	
	interp_barf_add_func(barf, opt_f);
	
	if(!barf->func)
		raiseErr("function to interpolate is missing. cannot add fvalues");
	
	SAFE_FREE(barf->f_values);
	SAFE_ALLOC(barf->f_values, barf->npoints);
	
	for(int i=0; i<barf->npoints; i++)
		barf->f_values[i] = barf->func(barf->points[i],barf->params);
}

static inline bool is_diff_too_small( double x1, double x2 )
{
	const double abs_diff = fabs( x1- x2 );
	const double tolerance = 10*DBL_EPSILON * ( 1 + fabs(x2) );
	return abs_diff <= tolerance;
}

double interp_barf_get_value( const BarFit *barf, double x )
{
	double num = 0;
	double den = 0;
	double temp = 0;
	const int nfitpoints = barf->npoints;
	
	for ( int i=0; i < nfitpoints; i++ )
	{
		if( is_diff_too_small( x, barf->points[i] ) )
			return barf->f_values[i];
		
		temp = barf->weights[i]/(x - barf->points[i]);
		num += temp * barf->f_values[i];
		den += temp;
	}
	
	return num/den;
}

void interp_barf_free(BarFit *barf)
{
	SAFE_FREE(barf->points);
	SAFE_FREE(barf->weights);
	SAFE_FREE(barf->f_values);
	SAFE_FREE(barf->params);
	barf->npoints=0;
	barf->func=NULL;
}

double* interp_cheb2_nodes_def( int num )
{
	double *nodes;
	SAFE_ALLOC(nodes,num);
	const double const_factor = PI/( num - 1 );
	for( int i = 0; i < num; i++ )
		nodes[i] = -cos( i  * const_factor );
	return nodes;
}

double* interp_cheb2_nodes( double x1, double x2, int num )
{
	double *nodes=interp_cheb2_nodes_def(num);

	const double half_length = ( x2 - x1 ) / 2;
	for( int i = 0; i < num; i ++ )
		nodes[i] = (nodes[i] + 1 ) * half_length + x1;
	return nodes;
}

double* interp_build_cheb2_weights( int num )
{
	double *weights;
	SAFE_ALLOC(weights,num);
	for( int i = 0; i < num; i ++ )
	{
		weights[i] = 1;
		if( i % 2 )
			weights[i] *= -1;
		if( i == 0 || i == num-1 )
			weights[ i ] /= 2;
	}
	return weights;
}

void interp_barf_init_cheb2_points( BarFit *barf, double x1, double x2, int npoints )
{
	barf->npoints = npoints;
	barf->points = interp_cheb2_nodes(x1,x2,npoints);
	barf->weights = interp_build_cheb2_weights( npoints );
}

BarFit interp_barf_new_cheb2_points(double x1, double x2, int npoints)
{
	BarFit barf={0};
	interp_barf_init_cheb2_points(&barf, x1, x2, npoints);
	return barf;
}
