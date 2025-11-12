#include "comp_physics.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

#define N_POINTS_INF_NORM 4000

// Calcola il numero di derivate continue di m

int m = 0;

double pown( double x, int n )
{
	double product = 1.0;
	for( int i = 0; i < n; i++ )
		product *=x;
	
	return product;
}

double f1(double x){ return pown( fabs( x ), m ); }

double inf_norm( const BarFit * barf )
{
	double max_diff = 0;
	const double length = 2;
	double curr_diff;
	double curr_point;
	for( int i = 0; i < N_POINTS_INF_NORM; i ++ )
	{
		curr_point = -1 + (length * i) /N_POINTS_INF_NORM;
		curr_diff = fabs( barf_get_value( barf, curr_point ) 
							- barf->func( curr_point ) );
		if( max_diff < curr_diff )
			max_diff = curr_diff;
	}
	
	
	return max_diff;
}

double one_test( Func_Ptr func, int n_nodes )
{
	BarFit barf = fit_cheb2_def( func, n_nodes );
	
	double sup_norm = inf_norm( &barf );
	
	free_vecD( &barf.points, &barf.f_values, &barf.weights );
	
	return sup_norm;
}
 
void test_cheb2( Func_Ptr func , const ArrayInt * n_nodes)
{
	VectorD sup_vec = init_vecD();
	for( int i = 0; i < n_nodes->length; i ++ )
		appendD( &sup_vec, one_test(func, n_nodes-> val[i]) );
	//printf( "Sup_norm is: %lf\n", sup_norm );
	
	VectorD n_nodesD = init_vecD();
	for( int i = 0; i < n_nodes->length; i++ )
		appendD( &n_nodesD, (double)n_nodes->val[i] );
		
	//printf( "%d %d %d\n", n_nodes->length, sup_vec.length, n_nodesD.length );
	
	plot_2vecs( &n_nodesD ,&sup_vec );
	eprint("ok2");
	free_vecD( &sup_vec, &n_nodesD );
}

int main()
{
	ArrayInt n_nodes = allocArrI( 10 );
	for( int i = 0; i < 10; i++ )
		n_nodes.val[i] = (i+1) * 10;
	
	//ArrayInt exp = allocArrI( 6 );
	for( int i = 0; i < 6; i++)
	{
		m = 1 + 2 * i;
		test_cheb2( f1, &n_nodes );
		eprint("ok1");
	}
	
	exit(EXIT_SUCCESS);
}

