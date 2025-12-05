//Es15_2.c

#include <stdio.h>
#include <stdlib.h>
#include <comp_physics.h>
#include <math.h>

#define DOMAIN_PLOT_SIZE 1000

double f1(double x, void * unused)
{
	UNUSED(unused);
	return 1/x/x - sin(x);
}

double f1_der(double x, void * unused)
{
	UNUSED(unused);
	return -2/x/x/x - cos(x);
}

Par_Func f1_p =     { f1,     NULL, 0 };
Par_Func f1_der_p = { f1_der, NULL, 0 };

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
	printf("Study of convergence for starting point %lf to root %.20lf:\n",
				debug->val[0], root);
	VectorD digits = init_vec_length( debug->length - 1 );
	for( int i = 0; i < digits.length; i ++ )
		digits.val[i] = -log10( fabs( debug->val[i] - root ) );
		
	VectorD ratios = init_vec_length( digits.length - 1 );
	for( int i = 0; i < ratios.length; i ++ )
		ratios.val[i] = digits.val[i + 1] / digits.val[i];
		
	std_print_vecD( &ratios );
}

void print_table( VectorD * starting_points, VectorD * roots )
{
	printf( "------------------------------------\n"
	        "----- Starting points and roots ----\n"
	        "------------------------------------\n"
	        "- starts ----------- roots ---------\n"
	        "------------------------------------\n");
	for( int i =0; i < starting_points->length; i++ )
		printf( "| %6g | %23.20g |\n", starting_points->val[i], roots->val[i] );
	printf( "------------------------------------\n");
}

void test( Par_Func f_p, Par_Func f_der_p, double x1, double x2, VectorD starting_points )
{
	VectorD debug = init_vecD();
	UNUSED(debug);

	const int nroots = starting_points.length;
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
	
	print_table( &starting_points, &roots );
	
	free_vecD( &debug, &roots );
}

int main()
{
	VectorD starting_points = init_vec_length(7);
	for( int i = 0; i < 7; i++ )
		starting_points.val[i] = (double)i + 1;
	test(f1_p, f1_der_p, 0.5, 10, starting_points);
	
	free_vecD(&starting_points);
	
	return 0;
}

