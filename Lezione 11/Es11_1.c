#include "comp_physics.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define PLOT_DOMAIN_LENGTH 100

double pow1(double x){return x;}
double pow2(double x){return x*x;}

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
	
	printf("Valori:\n");
	std_print_vecD( &fit.f_values );
	
	printf("\nPesi:\n");
	std_print_vecD(&fit.weights);
	printf("\n");
	
	// std_print_vecD( &domain );
	// std_print_vecD( &y_values );
	
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
	
	printf("Valori:\n");
	std_print_vecD( &fit.f_values );
	
	printf("\nPesi:\n");
	std_print_vecD(&fit.weights);
	printf("\n");
	
	// std_print_vecD( &domain );
	// std_print_vecD( &y_values );
	
	plot_2vecs( &domain, &y_values );
	
	free_vecD( &baric_fitter.points, &baric_fitter.weights,
			&fit.f_values, &fit.points, &fit.weights, &domain, &y_values);
}

int main()
{
	//func_test( pow1, 4, 1, 10 );
	func_test( pow2, 100, 0, 2 );
	func_test_eq(pow2, 100, 0 ,2);
	
	exit(EXIT_SUCCESS);
}
