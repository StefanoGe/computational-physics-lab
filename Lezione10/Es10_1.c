#include "comp_physics.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

double cosk( double x, double * k )
{
	return cos( x * *k );
}

Par_Func init_cosk(int k)
{
	Par_Func pf;
	pf.param_func_ptr = cosk;
	pf.params = malloc( sizeof( double ) );
	pf.params[0] = (double) k;
	return pf;
}

ArrayDouble generate_x_data( int n )
{
	return linspaceD( -1, 1, n + 1 );
}

ArrayDouble generate_y_data( ArrayDouble x_data, int k )
{
	ArrayDouble y_data = allocArrD( x_data.length );
	for( int i = 0; i < x_data.length; i++ )
		y_data.val[i] = cos( k * x_data.val[i] );
	return y_data;
}

ArrayDouble polynomial_interpolation( ArrayDouble x_data, ArrayDouble y_data)
{
	if( x_data.length != y_data.length )
		raiseErr( "x_data = %d and y_data = %d must have same size\n", 
			x_data.length, y_data.length );
	
	const int n_data = x_data.length;

	MatrixDouble A = allocMatD( n_data, n_data );
	
	// Set first col to 1 ( as it is the col of the constant function)
	for(int row = 0; row < n_data; row++ )
		A.val[row][0] = 1.0;
	
	for( int col = 1; col < n_data; col++ )
		for( int row = 0; row < n_data; row++ )
			A.val[ row ][ col ] = A.val[row][col - 1] * x_data.val[row];
			
	ArrayDouble result = solve_LU( A, y_data );
	
	freeMatD(A);
	
	return result;
}

void test( int n, int k )
{
	ArrayDouble x_data = generate_x_data(n);
	ArrayDouble y_data = generate_y_data( x_data, k );
	
	ArrayDouble coeffs = polynomial_interpolation(x_data, y_data);
	
	printf("For n = %d, k = %d\n", n, k);
	printArrDPar( coeffs, "%lf, " );
	putchar('\n');
	
	freeAllArrD( x_data, y_data, coeffs, NULL_ARR );
}



int main()
{
	for( int k = 10; k <= 100; k+=10 )
		for( int n = 4; n <= 40; n+=4 )
			test( n, k );
	
	exit(EXIT_SUCCESS);
}
