#include "comp_physics.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

// capire come valutare bene questa quantità

#define LOW_NUMBER 0.0001
#define DOMAIN_PLOT 1000

// Ricordati di studiare le proprietà della funzione
// e dimostrare che sia biunivoca (ha derivata maggiore di zero --> monotona)
// e ha come limite agli estremi di definizione +- inf

double f1(double z)
{
	return 1 / ( z*z -2*z + 2 );
}

VectorD cheb1_nodes_def(int num)
{
	const double common_factor = PI /(2 * num);
	VectorD nodes = init_vec_length(num);
	for(int i = 1; i <= num; i++)
	{
		nodes.val[i - 1] = -cos( (2 * i - 1) * common_factor );
	}

	
	//std_print_vecD( &nodes );
	
	return nodes;
}

VectorD cheb1_weights( int num )
{
	const double common_factor = PI / (num*2);
	VectorD weights = init_vec_length(num);
	printf("weights:\n");
	for(int i = 1; i <= num; i++)
	{
		weights.val[i - 1] = sin( (2 * i - 1) * common_factor );
		if( i % 2 == 1 )
			weights.val[i-1]*=(-1);
		printf("%lf ", weights.val[i-1]);
	}
	printf("\n");
	return weights;
}

VectorD cheb1_nodes( double x1, double x2, int num )
{
	VectorD nodes = cheb1_nodes_def( num );
	double half_length = ( x2 - x1 ) / 2;
	for( int i = 0; i < num; i ++ )
		nodes.val[i] = (nodes.val[i] + 1 ) * half_length + x1;
	return nodes;
}

BarFit fit_cheb1_def( Func_Ptr func, int n_nodes )
{
	BaricFitter fit_par;
	fit_par.points = cheb1_nodes_def( n_nodes );
	fit_par.weights = cheb1_weights( n_nodes );
	BarFit barf = bar_fir( &fit_par, func );
	free_vecD( &fit_par.points, &fit_par.weights );
	return barf;
}

BarFit fit_cheb1( Func_Ptr func, int n_nodes, double x1, double x2 )
{
	BaricFitter fit_par;
	fit_par.points = cheb1_nodes( x1, x2, n_nodes );
	fit_par.weights = cheb1_weights( n_nodes );
	BarFit barf = bar_fir( &fit_par, func );
	free_vecD( &fit_par.points, &fit_par.weights );
	return barf;
}

double unit_to_real_line( double x )
{
	// Here -1 < x < 1
	return 2 * x / ( 1 - x * x );
}

BarFit bar_fir_real( const BaricFitter * bar_fitter, Func_Ptr func )
{
	BarFit fit;
	fit.points = vec_cp( &bar_fitter->points );
	fit.weights = vec_cp( &bar_fitter->weights );
	fit.f_values = alloc_vecD( fit.points.size );
	fit.f_values.length = fit.points.length;
	for( int i = 0; i < fit.f_values.length; i ++ )
		fit.f_values.val[i] = func( unit_to_real_line( fit.points.val[i] ) );
	
	fit.func = func;
	
	return fit;
}

BarFit fit_cheb1_real( Func_Ptr func, int n_nodes )
{
	BaricFitter fit_par;
	fit_par.points = cheb1_weights(n_nodes);
	fit_par.weights = cheb1_weights( n_nodes );
	BarFit barf = bar_fir_real( &fit_par, func );
	free_vecD( &fit_par.points, &fit_par.weights );
	return barf;
}

//Ancora meglio: calcolare la funzione con una interpolazione!
//Ma aspe non sappiamo farla lol

double real_line_to_unit( double z )
{
	return ( fabs(z) < LOW_NUMBER ) ? 0 : ( sqrt( 1 + z * z ) -1 ) / z;
}

static inline bool is_diff_too_small( double x1, double x2 )
{
	const double abs_diff = fabs( x1- x2 );
	const double tolerance = DBL_EPSILON * ( 1 + fabs(x2) );
	return abs_diff <= tolerance;
}

double barf_get_value_real( const BarFit * barf, double z )
{
	double num = 0;
	double den = 0;
	const int length = barf->points.length;
	double temp = 0;
	
	double x = real_line_to_unit(z);
	
	
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

void test( int num )
{
	VectorD domain = vec_range( -6, 6, DOMAIN_PLOT );
	BaricFitter baric_fitter;
	baric_fitter.points = cheb1_nodes_def( num );
	baric_fitter.weights = cheb1_weights(num);
	
	BarFit barf = bar_fir_real( &baric_fitter, f1 );
	
	VectorD y_values = init_vec_length( DOMAIN_PLOT );
	VectorD true_values = init_vec_length(DOMAIN_PLOT);
	for( int i = 0; i < DOMAIN_PLOT; i++ )
	{
		y_values.val[i] = barf_get_value_real( &barf, domain.val[i] );
		true_values.val[i] = f1(domain.val[i]);
	}
	
	VectorD y_vecs [] = { y_values, true_values };
	
	plot_mult_vecs( &domain, y_vecs, 2, NULL_INFO );
	
	free_vecD( &baric_fitter.points, &baric_fitter.weights, &barf.f_values,
				&barf.points, &barf.weights, &y_values, &true_values);
}

void test_inv_func ( int num)
{
	VectorD cheb1 = cheb1_nodes_def( num );
	printf( "Cheb1 points:\n" );
	std_print_vecD(&cheb1);
	double curr;
	double unit;
	
	for( int i = 0; i < num; i++)
	{
		curr = unit_to_real_line( cheb1.val[i] );
		unit = real_line_to_unit(curr);
		printf( "%10lf   %10lf  %3lf\n", curr, unit, unit/cheb1.val[i] );
	}
	
	free_vecD( &cheb1 );
}

int main()
{
	//test_inv_func( 10 );
	
	
	test( 30 );
	
	
	exit(EXIT_SUCCESS);
}

