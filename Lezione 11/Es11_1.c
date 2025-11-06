#include "comp_physics.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

VectorD vec_range( double x1, double x2, int num )
{
	if( x2 < x1 )
		raiseErr( "given x2 < x1" );
	const double spacing = (x2- x1) / (num - 1);
	VectorD vec = init_vecD();
	for( int i = 0; i < num; i++ )
		appendD( &vec, i * spacing + x1 );
	
	return vec;
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
	
	std_print_vecD( &fit.f_values );
	
	free_vecD( &points, &baric_fitter.points, &baric_fitter.weights,
			&fit.f_values, &fit.points, &fit.weights);
}

int main()
{
	func_test( log, 3, 1, 10 );
	
	VectorD domain = vec_range( 0, 1, 1000 );
	
	eprint("ok");
	
	plot_func( &domain, log );
	
	
	exit(EXIT_SUCCESS);
}
