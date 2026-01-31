#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "comp_physics.h"

/*
static inline VectorD arrD_as_VecD( ArrayDouble arr )
{
	VectorD vec;
	vec.length = arr.length;
	vec.val = arr.val;
	vec.size = vecD_necessary_size( vec.length );
	vec.val = realloc( vec.val, vec.size * sizeof(double) );
	return vec;
}
*/

static inline VectorD build_weights( const VectorD * points )
{
	const int length = points->length;
	VectorD weights = init_vec_length( length );
	double prod = 1;
	
	for( int w = 0; w < length; w++ )
	{
		prod = 1;
		for( int j = 0; j < length; j++ )
			if( w != j )
				prod *= ( points->val[w] - points -> val[j] );
				
		weights.val[w] = 1 / prod;
	}
	
	return weights;
}

BaricFitter init_bar_fitter( const VectorD * points )
{
	BaricFitter barf;
	barf.points = vec_cp( points );
	barf.weights = build_weights( points );
	return barf;
}

static inline VectorD eq_weights( int n )
{
	VectorD weights = init_vec_length( n );
	weights.val[ 0 ] = 1;
	for( int i = 1; i < n; i++ )
		weights.val[i] = weights.val[ i -1 ] * ( i - 1 - n ) / i;
	return weights;
}

BaricFitter bar_fitter_eq_init( double x1, double x2, int n_points )
{
	BaricFitter barf;
	barf.points = vec_range( x1, x2, n_points );
	barf.weights = eq_weights(n_points);
	return barf;
}

BarFit bar_fir( const BaricFitter * bar_fitter, Func_Ptr func )
{
	BarFit fit;
	fit.points = vec_cp( &bar_fitter->points );
	fit.weights = vec_cp( &bar_fitter->weights );
	fit.f_values = alloc_vecD( fit.points.size );
	fit.f_values.length = fit.points.length;
	for( int i = 0; i < fit.f_values.length; i ++ )
		fit.f_values.val[i] = func( fit.points.val[i] );
	
	fit.func = func;
	
	return fit;
}

static inline bool is_diff_too_small( double x1, double x2 )
{
	const double abs_diff = fabs( x1- x2 );
	const double tolerance = DBL_EPSILON * ( 1 + fabs(x2) );
	return abs_diff <= tolerance;
}

double barf_get_value( const BarFit * barf, double x )
{
	double num = 0;
	double den = 0;
	const int length = barf->points.length;
	
	double temp = 0;
	
	for ( int i=0; i < length; i++ )
	{
		if( is_diff_too_small( x, barf -> points.val[i] ) )
			return barf->f_values.val[i];
		
		temp = barf->weights.val[i]/( x - barf -> points.val[i] );
		num += temp * barf->f_values.val[i];
		den += temp;
	}
	
	return num/den;
}

VectorD cheb2_nodes_def( int num )
{
	VectorD nodes = init_vec_length( num );
	const double const_factor = PI/( num - 1 );
	for( int i = 1; i <= num; i++ )
		nodes.val[i-1] = -cos( ( i - 1 ) * const_factor );
	return nodes;
}

VectorD cheb2_nodes( double x1, double x2, int num )
{
	VectorD nodes = cheb2_nodes_def( num );
	double half_length = ( x2 - x1 ) / 2;
	for( int i = 0; i < num; i ++ )
		nodes.val[i] = (nodes.val[i] + 1 ) * half_length + x1;
	return nodes;
}

VectorD cheb2_weights( int num )
{
	VectorD nodes = init_vec_length( num );
	for( int i = 1; i <= num; i ++ )
	{
		nodes.val[i - 1] = 1;
		if( i % 2 )
			nodes.val[i - 1] *= (-1);
		if( i == 1 || i == num )
			nodes.val[ i - 1 ] /= 2;
	}
	return nodes;
}

BarFit fit_cheb2_def( Func_Ptr func, int n_nodes )
{
	BaricFitter fit_par;
	fit_par.points = cheb2_nodes_def( n_nodes );
	fit_par.weights = cheb2_weights( n_nodes );
	BarFit barf = bar_fir( &fit_par, func );
	free_vecD( &fit_par.points, &fit_par.weights );
	return barf;
}

BarFit fit_cheb2( Func_Ptr func, int n_nodes, double x1, double x2 )
{
	BaricFitter fit_par;
	fit_par.points = cheb2_nodes( x1, x2, n_nodes );
	fit_par.weights = cheb2_weights( n_nodes );
	BarFit barf = bar_fir( &fit_par, func );
	free_vecD( &fit_par.points, &fit_par.weights );
	return barf;
}



