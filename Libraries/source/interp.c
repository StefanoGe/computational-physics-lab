#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "comp_physics.h"

/*
typedef struct BaricFit{
	VectorD points;
	VectorD weights;
	VectorD coeffs;
	
};
*/

static inline double get_kth_vector( const VectorD * points, int step )
{
	double prod = 1;
	for( int i = 0; i < step; i++ )
		prod *= points->val[step] - points->val[i];
		
	return prod;
}

VectorD build_weights( const VectorD * points )
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
	
	VectorD weight_vec;
	weight_vec.length = length;
	int size = DEFAULT_VEC_SIZE;
	while( (size *=2) < length );
	weight_vec.size = size;
	weights = realloc( weights, size * sizeof(double) );
	return weight_vec;
}

