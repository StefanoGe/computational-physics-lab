#include <stdio.h>
#include "comp_physics.h"
#include <math.h>

#define DOMAIN_PLOT_SIZE 1000

double f1( double x, void * unused )
{
	UNUSED(unused);
	return x + 1;
}

void test1()
{
	Par_Func fnc_obj;
	fnc_obj.param_func_ptr = f1;
	fnc_obj.params = NULL;
	fnc_obj.nparams = 0;
	
	printf("%lf", root_bis( fnc_obj, -2.5, 0, 1e-10 , NULL) );
}

double f2(double x, void * unused)
{
	UNUSED(unused);
	return x*x*x - 7*x*x + 14*x -6;
}

void plot_graph( double root )
{
	VectorD domain = vec_range( 0, 1, DOMAIN_PLOT_SIZE );
	VectorD y_values = init_vec_length( DOMAIN_PLOT_SIZE );
	for( int i =0; i < DOMAIN_PLOT_SIZE; i++ )
		y_values.val[i] = f2( domain.val[i] , NULL );
	
	FILE * gp = gp_open();
	gp_term_def( gp, "Graph of particular cubic function" );
	fprintf(gp, "set grid nopolar\n");
	gp_axes_labels( gp, "x", "f(x)" );
	char * style2 = "points";
	char * styles[2];
	styles[0] = 0;
	styles[1] = style2;
	gp_set_plot( gp, 2, NULL, styles );
	gp_prt_carr( gp, domain.val, y_values.val, domain.length );
	fprintf( gp, "%g %g\ne\n", root, 0.0 );
	gp_end(gp);
	
	free_vecD( &domain, &y_values );
}

void plot_graph_conv(VectorD * log_abs_conv)
{
	FILE * gp = gp_open();
	gp_term_def( gp, "Conergence of bisection method" );
	fprintf(gp, "set grid nopolar\n");
	gp_axes_labels( gp, "n", "-log..." );
	char * styles[] = {"lp"};
	gp_set_plot( gp, 1, NULL, styles );
	VectorD domain = vec_range( 1, log_abs_conv->length, log_abs_conv->length );
	gp_prt_carr( gp, domain.val, log_abs_conv->val, log_abs_conv->length );
	gp_end(gp);
}

void test2(double x1, double x2)
{
	Par_Func fnc_obj;
	fnc_obj.nparams = 0;
	fnc_obj.param_func_ptr = f2;
	fnc_obj.params = NULL;
	
	VectorD debug = init_vecD();
	
	double root = root_bis( fnc_obj, x1, x2, 0, &debug );
	
	printf( "Found root %.20lf\n", root );
	
	print_vecD( &debug, "%.20lf\n", stdout );
	
	VectorD abs_err = init_vecD();
	VectorD log_abs_error = init_vecD();
	
	for( int i = 0; i < debug.length; i++ )
	{
		appendD( &abs_err, fabs( debug.val[i] - root ) );
		appendD( &log_abs_error, fabs( -log10( abs_err.val[i] ) ) );
	}
	
	//plot_graph( root );
	plot_graph_conv( &log_abs_error );
	
	free_vecD( &debug );
}

int main()
{
	test2(0,1);
	test2(0.4, 1);
	test2(0.25, 0.66);
	
	exit(EXIT_SUCCESS);
}

// ci dovrebbe uscire 0.5 perche' ogni volta dividiamo l'intervallo in 2
