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
 
VectorD test_cheb2( Func_Ptr func , const ArrayInt * n_nodes)
{
	VectorD sup_vec = init_vecD();
	for( int i = 0; i < n_nodes->length; i ++ )
		appendD( &sup_vec, one_test(func, n_nodes-> val[i]) );
	
	return sup_vec;
}

int main()
{
	ArrayInt n_nodes = allocArrI( 10 );
	for( int i = 0; i < 10; i++ )
		n_nodes.val[i] = (i+1) * 10;
		
	VectorD n_nodesD = init_vecD();
	for( int i = 0; i < n_nodes.length; i++ )
		appendD( &n_nodesD, (double)n_nodes.val[i] );
	
	VectorD sup_norms [6];
	char * labels[6];
	
	//ArrayInt exp = allocArrI( 6 );
	for( int i = 0; i < 6; i++)
	{
		m = 1 + 2 * i;
		labels[i] = malloc( 100 * sizeof(char) );
		sprintf( labels[i], "m = %d", m );
		sup_norms[i] = test_cheb2( f1, &n_nodes );
	}
	
	FILE * gp = gp_open();
	
	char * title = malloc( sizeof(char) * 100 );
	sprintf(title, "Plot finale");
	
	fprintf(gp, "set logscale y\n");
	fprintf(gp, "set key box\n");
	
	gp_term_def( gp, title );
	gp_axes_labels( gp, "n", "inf-norm error" );
	gp_set_plot(gp, 6, labels);
	for( int i = 0; i < 6; i++ )
	{
		std_print_vecD( sup_norms+i );
		gp_prt_carr( gp, n_nodesD.val, sup_norms[i].val, n_nodesD.length );
	}
	
	gp_end(gp);
	free(title);
	
	exit(EXIT_SUCCESS);
}



