// Test of trapezoidal rule of integration and study of its convergence

#include <stdio.h>
#include <math.h>
#include "comp_physics.h"

#define MAX_NODE_POWER 10
#define DO_ALL false

double f1( double x, void * unused )
{
	UNUSED(unused);
	return x * log(1 + x);
}

Par_Func f1_p = { f1, NULL, 0 };

double f2( double x, void * unused )
{
	UNUSED(unused);
	return x*x*atan(x);
}

Par_Func f2_p = { f2, NULL, 0 };

double f3( double x, void * unused )
{
	UNUSED(unused);
	return exp(x) * cos(x);
}

Par_Func f3_p = { f3, NULL, 0 };

double f4( double x, void * unused )
{
	UNUSED(unused);
	const double A = sqrt( 2 + x * x );
	return atan( A ) / ( 1 + x * x )/A;
}

Par_Func f4_p = { f4, NULL, 0 };

double f5( double x, void * unused )
{
	UNUSED(unused);
	if (x == 0)
		return 0;
	return sqrt(x) * log(x);
}

Par_Func f5_p = { f5, NULL, 0 };

double f6( double x, void * unused )
{
	UNUSED(unused);
	return sqrt( 1 - x * x );
}

Par_Func f6_p = { f6, NULL, 0 };

int pow_int( int base, int exp )
{
	int prod = 1;
	for( int i = 0; i < exp; i++ )
		prod*=base;
	return prod;
}

double abs_err( double value, double true_value )
{
	return fabs(value - true_value);
}

void plot_graph_conv(const VectorD * n_node_list, const VectorD * errors )
{
	FILE * gp = gp_open();
	gp_term_def( gp, "Convergence of trapezoidal rule" );
	fprintf(gp, "set grid nopolar\n");

	fprintf(gp, "set logscale x\n");	
	fprintf(gp, "set logscale y\n");

	gp_axes_labels( gp, "n nodes", "Error" );
	char * styles[] = {"lp"};

	gp_set_plot( gp, 1, NULL, styles );
	gp_prt_carr( gp, n_node_list->val, errors->val, n_node_list->length );
	gp_end(gp);
}

void integration_study( const Par_Func * fnc, double x1, double x2, double true_value )
{
	VectorD errors = init_vecD();
	VectorD n_node_list = init_vecD();
	printf("True value: %g\n", true_value);
	for( int i = 1; i <= MAX_NODE_POWER; i++ )
	{
		int n_nodes = 10 * pow_int(2, i);
		double integral = int_trap( fnc, x1, x2, n_nodes );
		double error = true_value - integral;
		//Absolute error or no? Yes for the log plot
		appendD( &errors, fabs(error) );
		appendD( &n_node_list, (double) n_nodes );
		printf( "n_nodes: %6d. Error: %e.\n", n_nodes, error );
	}
	eprint("ok before tmp");
	tmplot_2vecs( "trap_rule_plot", &n_node_list, &errors );
	
	free_vecD(&errors, &n_node_list);
	// Aggiungere interpolazione
}

int main()
{
	integration_study( &f1_p, 0, 1, 0.25 );
	if(DO_ALL)
	{
		integration_study( &f2_p, 0, 1, (PI - 2 + 2 * log(2)) / 12 );	
		integration_study( &f3_p, 0, PI/2, ( exp( PI/2 ) - 1 ) / 2 );
		integration_study( &f4_p, 0, 1, 5 * PI*PI/96 );
		integration_study( &f5_p, 0, 1, -4.0/9 );
		integration_study( &f6_p, 0, 1, PI/4 );
	}
	
	
	exit(EXIT_SUCCESS);
}

