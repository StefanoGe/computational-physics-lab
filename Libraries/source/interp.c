#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "comp_physics.h"



static inline double get_kth_vector( const VectorD * points, int step )
{
	double prod = 1;
	for( int i = 0; i < step; i++ )
		prod *= points->val[step] - points->val[i];
		
	return prod;
}

static inline VectorD cp_carr_to_vec( double * carr, int length )
{
	VectorD vec;
	
	vec.length = length;
	int size = DEFAULT_VEC_SIZE;
	while( (size *=2) < length );
	vec.size = size;
	carr = realloc( carr, size * sizeof(double) );
	vec.val = carr;
	
	return vec;
}

static inline VectorD build_weights( const VectorD * points )
{
	const int length = points->length;
	double * weights = malloc( length * sizeof(double) );
	if( !weights )
		raiseErr( "Failed to allocate memory!" );
	
	for( int step = 0; step < length; step ++ )
	{
		weights[ step ] = get_kth_vector( points, step );
		for( int w = 0; w < step; w++ )
		{
			if( w != step )
				weights[ w ] *= weights[w]*(points->val[w] - points->val[step]);
		}
	}
	
	for( int i =0; i < length; i ++ )
		weights[i] = 1/weights[i];
	
	return cp_carr_to_vec( weights, length );
}

BaricFitter init_bar_fitter( const VectorD * points )
{
	BaricFitter barf;
	barf.points = vec_cp( points );
	barf.weights = build_weights( points );
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
	
	return fit;
}


