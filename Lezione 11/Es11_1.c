#include "comp_physics.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define PLOT_DOMAIN_LENGTH 100

double pow1(double x){return x;}
double pow2(double x){return x*x;}


VectorD fitted_y_values( VectorD * domain, Func_Ptr func, 
							int n_nodes, double x1, double x2 )
{
	VectorD points = vec_range( x1, x2, n_nodes );
	BaricFitter baric_fitter = init_bar_fitter(&points);
	BarFit fit = bar_fir( &baric_fitter, func );

	VectorD y_values = init_vec_length( PLOT_DOMAIN_LENGTH );
	
	for(int i = 0; i < PLOT_DOMAIN_LENGTH; i++)
		y_values.val[i] = barf_get_value( &fit, domain->val[i] );
		
	free_vecD(&baric_fitter.points, &baric_fitter.weights, &fit.f_values, &fit.points,
	&fit.weights);
	return y_values;
}

void func_test_eq_2( Func_Ptr func, int n_nodes[], int n_tests, double x1, double x2 )
{
	VectorD domain = vec_range( x1, x2, PLOT_DOMAIN_LENGTH );
	
	VectorD * y_values = malloc( (n_tests + 1) * sizeof(VectorD) );
	for(int i = 0; i < n_tests; i++)
		y_values[i] = fitted_y_values( &domain, func, n_nodes[i],
							x1, x2);
	
	VectorD true_values = init_vec_length( PLOT_DOMAIN_LENGTH );
	for( int i = 0; i < true_values.length; i++ )
		true_values.val[i] = func( domain.val[i] );
	y_values[n_tests] = true_values;
	
	PlotInfo plot_info = plot_info_init( n_tests );
	for( int i = 0; i < n_tests; i++ )
	{
		plot_info.labels[i] = malloc( sizeof(char) * 20 );
		sprintf( plot_info.labels[i], "n = %d", n_nodes[i] );
	}
	plot_info.labels[n_tests] = "True values";
	
	plot_mult_vecs( &domain, y_values, n_tests + 1, plot_info );
	
	for(int i = 0; i < n_tests; i++)
		free_vecD( &y_values[i] );
	free_vecD( &domain);
}


int main()
{
	int n_nodes_1 []= {2, 3, 4};
	func_test_eq_2( log, n_nodes_1, 3, 1, 10 );
	func_test_eq_2( tanh, n_nodes_1, 3, -3, 2 );
	func_test_eq_2( cosh, n_nodes_1, 3, -1, 3 );
	int n_nodes_2[] = {3, 5, 7};
	func_test_eq_2( fabs, n_nodes_2, 3, -2, 1 );
	
	exit(EXIT_SUCCESS);
}

void first_test()
{
	VectorD points = build_vecD( 4, 0.0, 1.0, 2.0, 3.0 );
	
	std_print_vecD( &points );
	
	free_vecD(&points );
}

void func_test( Func_Ptr func, int n_nodes, double x1, double x2 )
{
	VectorD points = vec_range( x1, x2, n_nodes );
	BaricFitter baric_fitter = init_bar_fitter(&points);
	BarFit fit = bar_fir( &baric_fitter, func );
	
	VectorD domain = vec_range( x1, x2, PLOT_DOMAIN_LENGTH );
	VectorD y_values = init_vec_length( PLOT_DOMAIN_LENGTH );
	for(int i = 0; i < PLOT_DOMAIN_LENGTH; i++)
		y_values.val[i] = barf_get_value( &fit, domain.val[i] );
	
	plot_2vecs( &domain, &y_values );
	
	free_vecD( &points, &baric_fitter.points, &baric_fitter.weights,
			&fit.f_values, &fit.points, &fit.weights, &domain, &y_values);
}

void func_test_eq( Func_Ptr func, int n_nodes, double x1, double x2 )
{
	BaricFitter baric_fitter = bar_fitter_eq_init(x1, x2, n_nodes );
	BarFit fit = bar_fir( &baric_fitter, func );
	
	VectorD domain = vec_range( x1, x2, PLOT_DOMAIN_LENGTH );
	VectorD y_values = init_vec_length( PLOT_DOMAIN_LENGTH );
	for(int i = 0; i < PLOT_DOMAIN_LENGTH; i++)
		y_values.val[i] = barf_get_value( &fit, domain.val[i] );
	
	plot_2vecs( &domain, &y_values );
	
	free_vecD( &baric_fitter.points, &baric_fitter.weights,
			&fit.f_values, &fit.points, &fit.weights, &domain, &y_values);
}
