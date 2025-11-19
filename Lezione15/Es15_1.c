#include <stdio.h>
#include <stdlib.h>
#include <comp_physics.h>
#include <math.h>

#define DOMAIN_PLOT_SIZE 1000

double f1( double x, double * unused )
{
	UNUSED(unused);
	return x*x - exp( -x );
}

double f1_der( double x, double * unused )
{
	UNUSED(unused);
	return 2 * x + exp(-x);
}

Par_Func f1_p =     { f1,     NULL, 0 };
Par_Func f1_der_p = { f1_der, NULL, 0 };

double f2( double x, double * unused )
{
	UNUSED(unused);
	return tan(x) - 2 * x;
}

double f2_der( double x, double * unused )
{
	UNUSED(unused);
	return 1/cos(x)/cos(x) - 2 ;
}

Par_Func f2_p =     { f2,     NULL, 0 };
Par_Func f2_der_p = { f2_der, NULL, 0 };

double f3( double x, double * unused )
{
	UNUSED(unused);
	return exp(x + 1) - x - 2;
}

double f3_der( double x, double * unused )
{
	UNUSED(unused);
	return exp(x+1) - 1;
}

Par_Func f3_p =     { f3,     NULL, 0 };
Par_Func f3_der_p = { f3_der, NULL, 0 };

void plot_graph( Par_Func f_p, double x1, double x2, const VectorD * roots)
{
	const int nroots = roots->length;
	VectorD domain = vec_range( x1, x2, DOMAIN_PLOT_SIZE );
	VectorD y_values = vec_par_func( &domain, &f_p );
	VectorD zeroes = vec_full( nroots, 0 );
	
	
	FILE * gp = gp_open();
	gp_term_def( gp, "Graph and roots of " );
	fprintf(gp, "set grid nopolar\n");
	gp_axes_labels( gp, "x", "f(x)" );
	char * style2 = "points";
	char * styles[2];
	styles[0] = NULL;
	styles[1] = style2;
	gp_set_plot( gp, 2, NULL, styles );
	gp_prt_carr( gp, domain.val, y_values.val, domain.length );
	gp_prt_carr( gp, roots->val, zeroes.val, nroots );
	gp_end(gp);
	
	free_vecD( &domain, &y_values, &zeroes );
}

void convergence_study(double root, const VectorD * debug)
{	
	VectorD digits = init_vec_length( debug->length - 1 );
	for( int i = 0; i < digits.length; i ++ )
		digits.val[i] = -log10( fabs( debug->val[i] - root ) );
		
	VectorD ratios = init_vec_length( digits.length - 1 );
	for( int i = 0; i < ratios.length; i ++ )
		ratios.val[i] = digits.val[i + 1] / digits.val[i];
		
	std_print_vecD( &ratios );
}

//NOTE: DETERMINE THE C CONSTANT!

void test( Par_Func f_p, Par_Func f_der_p, double x1, double x2, int nroots, ... )
{
	VectorD debug = init_vecD();
	UNUSED(debug);
	
	va_list args;
	va_start( args, nroots );
	
	VectorD starting_points = init_vec_length(nroots);
	for(int i =0; i < nroots; i++)
		starting_points.val[i] = va_arg( args, double );

	va_end(args);

	VectorD roots = init_vec_length(nroots);
	
	for( int i = 0; i < nroots; i++)
	{
		roots.val[i] = root_newt( &f_p, &f_der_p, starting_points.val[i],
					DEF_TOL, DEF_TOL, &debug );
		printf("Intermediate values for finding the root %lf of index %d:\n", 
				roots.val[i], i);
		print_vecD( &debug, "%.20lf\n", stdout );
		convergence_study( roots.val[i], &debug );
	}
	
	plot_graph( f_p, x1, x2, &roots );
	
	free_vecD( &debug, &roots, &starting_points );
}

int main()
{
	test( f1_p, f1_der_p, -2, 2, 1, 0.0 );
	test( f2_p, f2_der_p, -0.2, 1.4, 2, -0.1, 1.0 ); 
	//oss: convergenza del primo zero e' rapidissima poiche' la funzione e'
	//molto lineare
	test( f3_p, f3_der_p, -2, 2, 1, 0.0 );
	//Qui la convergenza e' molto lenta poich' la funzione ha derivata
	//0 nel punto cercato
	
	exit(EXIT_SUCCESS);
}

/* OLD
 * 
 * void test1( double x1, double x2  )
{
	VectorD debug = init_vecD();
	UNUSED(debug);
	
	double root = root_newt( &f1_p, &f1_der_p, 0, DEF_TOL, DEF_TOL, &debug );
	print_vecD( &debug, "%.20lf\n", stdout );
	
	VectorD roots = init_vec_length( 1 );
	roots.val[0] = root;
	
	plot_graph( f1_p, x1, x2, &roots );
	
	free_vecD( &debug, &roots );
}

void test2( double x1, double x2  )
{
	VectorD debug = init_vecD();
	UNUSED(debug);
	
	VectorD starting_points = build_vecD( 2, 0.0, 1.0 );
	const int nroots = starting_points.length;
	VectorD roots = init_vec_length(nroots);
	
	for( int i = 0; i < nroots; i++)
	{
		roots.val[i] = root_newt( &f2_p, &f2_der_p, starting_points.val[i],
					DEF_TOL, DEF_TOL, &debug );
		printf("Intermediate values for finding the root %lf of index %d:\n", 
				roots.val[i], i);
		print_vecD( &debug, "%.20lf\n", stdout );
	} 
	
	plot_graph( f2_p, x1, x2, &roots );
	
	free_vecD( &debug, &roots, &starting_points );
}
* 
*/


