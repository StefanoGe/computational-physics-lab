#include "comp_physics.h"
#include <stdio.h>
#include <math.h>

#define N_XS 100
#define NORMALIZATION_CONSTANT 2006.787453080206


double h( double t )
{
	return exp( sin( 4 * t ) );
}

void test()
{
	ArrayDouble x_data = linspaceD( 0,1, N_XS );
	ArrayDouble y_data = build_from( x_data, h );
	
	ArrayDouble coeffs = linear_ls_pol_fitting( x_data, y_data, 14 );
	
	printArrDPar( coeffs, "%.16lf " );
	
	free_arrD( x_data, y_data, coeffs );
}

int main()
{
	test();
	
	
	exit(EXIT_SUCCESS);
}
