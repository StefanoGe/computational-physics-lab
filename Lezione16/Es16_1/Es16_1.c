// Es16_1.c 
// Secant method implementation and study of convergence

#include <stdio.h>
#include <stdlib.h>
#include <comp_physics.h>
#include <math.h>

#define DOMAIN_PLOT_SIZE 1000

double f1( double x, void * unused )
{
	UNUSED(unused);
	return x*x - exp( -x );
}

Par_Func f1_p = { f1,     NULL, 0 };

double f2( double x, void * unused )
{
	UNUSED(unused);
	return tan(x) - 2 * x;
}

Par_Func f2_p = { f2,     NULL, 0 };

double f3( double x, void * unused )
{
	UNUSED(unused);
	return exp(x + 1) - x - 2;
}

Par_Func f3_p = { f3,     NULL, 0 };

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

void test( Par_Func f_p, double x1, double x2,
					double xtol, double ftol, int nroots, ... )
{
	VectorD debug = init_vecD();
	UNUSED(debug);
	
	va_list args;
	va_start( args, nroots );
	
	VectorD starting_points1 = init_vec_length(nroots);
	VectorD starting_points2 = init_vec_length(nroots);
	
	for(int i =0; i < nroots; i++)
	{
		starting_points1.val[i] = va_arg( args, double );
		starting_points2.val[i] = va_arg( args, double );		
	}

	va_end(args);

	VectorD roots = init_vec_length(nroots);
	
	for( int i = 0; i < nroots; i++)
	{
		int exit_status = -2;
		roots.val[i] = root_sec( &f_p, starting_points1.val[i], 
							starting_points2.val[i], xtol, ftol, &debug,
							&exit_status);
		
		if(exit_status == -1)
		{
			eprint("root_finding failed. debug:\n");
			print_vecD( &debug, "%.20lf, ", stderr );
		}
		
		if( exit_status != -2 && exit_status != -1 )
		{
			printf("Intermediate values for finding the root %lf of index %d:\n", 
					roots.val[i], i);
			print_vecD( &debug, "%.20lf\n", stdout );
			convergence_study( roots.val[i], &debug );
		}
		

	}
	 
	plot_graph( f_p, x1, x2, &roots );
	
	free_vecD( &debug, &roots, &starting_points1, &starting_points2 );
}

int main()
{
	test( f1_p, -2, 2, DEF_TOL, DEF_TOL, 1, -2.0, 2.0 );
	test( f2_p, -0.2, 1.4, DEF_TOL, DEF_TOL, 2, -0.2, 0.5, 0.9, 1.4 );
	test( f3_p, -2, 2, DEF_TOL, DEF_TOL, 1,-2.0, 2.0 );
	
	exit(EXIT_SUCCESS);
}

