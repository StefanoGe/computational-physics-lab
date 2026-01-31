#include "comp_physics.h"
#include <math.h>

#define MAX_DEGREE_TEST 80

double leg2roots[2];
double leg3roots[3];
double leg4roots[4];
double leg5roots[5];
double *leg_roots[]= {leg2roots, leg3roots, leg4roots, leg5roots};

void setup_leg_roots(){
	leg2roots[1]=-sqrt(3)/3;
	leg2roots[0]= sqrt(3)/3;
	leg3roots[0]= sqrt(15)/5;
	leg3roots[2]=-sqrt(15)/5;
	leg3roots[1]=0;
	leg4roots[3]=-sqrt( 525 + 70 * sqrt(30) )/35;
	leg4roots[0]= sqrt( 525 + 70 * sqrt(30) )/35;
	leg4roots[2]=-sqrt( 525 - 70 * sqrt(30) )/35;
	leg4roots[1]= sqrt( 525 - 70 * sqrt(30) )/35;

}

void test_roots(){
	
	VectorD debug;
	
	for( int deg = 1; deg <= MAX_DEGREE_TEST; deg++ )
	{
		printf(" n = %d\n", deg);
		for( int root = 1; root <= deg; root++ )
		{
			double curr_root = legendre_root(deg, root, &debug);
			printf( " k = %2d: %.16g\n", root, curr_root );
			if( 2 <=deg && deg <=4 )
				printf(" diff = %16g\n", fabs(curr_root-leg_roots[deg-2][root-1]) );
			printf("steps to find the root:\n");
			for( int i = 0; i < 100 && i < debug.length; i++) 
				printf( "%.20g  ", debug.val[i] );
			printf("\n");
		}
		printf(" -------------------------------------\n");
	}
	
	free_vecD(&debug);
}

void plot_polynomials()
{
	VectorD xaxis = vec_range( -1, 1, 100);
	MatrixDouble yvalues = allocMatD( 2, 100 );
	for( int i = 0; i < 100; i++)
	{
		double x = xaxis.val[i];
		yvalues.val[0][i] = legendre( x, 2 );
		yvalues.val[1][i] = ( 3*x*x - 1 ) / 2;
	}
	
	tmplot_carrs( "plot_pol", xaxis.val, yvalues.val, 100, 2 );
	
	freeMatD(yvalues);
	free_vecD(&xaxis);
}


void plot_polynomials24()
{
	int n_points = 1000;
	
	VectorD xaxis = vec_range( -1, 1, n_points);
	MatrixDouble yvalues = allocMatD( 1, n_points );
	for( int i = 0; i < n_points; i++)
	{
		double x = xaxis.val[i];
		yvalues.val[0][i] = legendre( x, 24 );
	}
	
	tmplot_carrs( "plot_pol", xaxis.val, yvalues.val, n_points, 1 );
	
	freeMatD(yvalues);
	free_vecD(&xaxis);
}

// Remark: it was necessary to set 1000 DBL epsilon in ftol. A safer 
// implementation would set them higher. 
// Remark: legendre memoization

int main()
{
	setup_leg_roots();
//	eprint("setup ok");
//	plot_polynomials24();
//	eprint("plot ok");
	test_roots();

	
	exit(EXIT_SUCCESS);
}

