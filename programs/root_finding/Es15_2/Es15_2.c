//Es15_2.c
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rootf.h"
#include "explot.h"

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

void arr_line_between_points(Array *x, Array *y, double x1, double y1,
			double x2, double y2, int n)
{
	arr_init(x,n);
	arr_init(y,n);
	const double factor_x = (x2-x1)/(n-1);
	const double factor_y = (y2-y1)/(n-1);
	for(int i=0; i<n; i++)
	{
		ARRP(x,i) = x1 + i*factor_x;
		ARRP(y,i) = y1 + i*factor_y;
	}
}

/*
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
*/

void plot_graph( ParamFuncPtr f, double x1, double x2, const Array *roots,
	const Array *initial_guesses)
{
	
	Array domain = arr_linspace(x1,x2,DOMAIN_PLOT_SIZE);
	Array y_values = arr_map_par(&domain,f,NULL);
	
	SeriesSpec graph = eplot_init_series(domain.data,y_values.data,domain.size,
					"l","f(x)",NULL);
	
	Array y_roots = arr_map_par(roots,f,NULL);
	
	SeriesSpec root_rep = eplot_init_series(roots->data,y_roots.data,roots->size,
					"p pt 2","Roots",NULL);
					
	Array initial_guesses_y = arr_new(initial_guesses->size);
	arr_setv(&initial_guesses_y,-1);
	
	SeriesSpec init_guesses_series = eplot_init_series(initial_guesses->data,
		initial_guesses_y.data,initial_guesses->size,"p pt 4","Initial guesses",
		NULL);
	
	Array *x_lines;
	SAFE_ALLOC(x_lines,roots->size);
	Array *y_lines;
	SAFE_ALLOC(y_lines,roots->size);
	
	SeriesSpec *series;
	SAFE_ALLOC(series,3+roots->size);
	
	for(int i=0; i<roots->size; i++)
	{
		arr_line_between_points(x_lines+i,y_lines+i,
			ARRP(initial_guesses,i),ARR(initial_guesses_y,i),
			ARRP(roots,i),0,100);
		
		series[i+3] = eplot_init_series(x_lines[i].data,
			y_lines[i].data,100,"l lw 1",NULL,
			NULL);
	}
	
	Global gb_settings = eplot_init_global("Plot of function x^{-2} - sin(x)",
						"x","y",NULL,"newt_root_discussion",NULL,NULL,
						NULL);
	
	series[0] = graph;
	series[1] = root_rep;
	series[2] = init_guesses_series;
	
	eplot_multi(series,3+roots->size,&gb_settings);
	
	arr_free(&domain);
	arr_free(&y_values);
	arr_free(&y_roots);
	arr_free(&initial_guesses_y);
	for(int i=0; i<roots->size; i++)
	{
		arr_free(x_lines+i);
		arr_free(y_lines+i);
	}
	SAFE_FREE(series);
	SAFE_FREE(x_lines);
	SAFE_FREE(y_lines);
}

void test( ParamFuncPtr f, ParamFuncPtr der, 
			double x1, double x2, const Array *starting_points )
{
	Array debug = {0};

	UNUSED(x1);
	UNUSED(x2);

	const int nroots_searched = starting_points->size;
	Array roots = arr_new(nroots_searched);
	
	for( int i = 0; i < nroots_searched; i++)
	{
		ARR(roots,i) = root_newt( f, NULL,der,NULL,ARRP(starting_points,i),
					DEF_TOL, DEF_TOL, &debug );
		printf("Intermediate values for finding the root %lf of index %d:\n", 
				ARR(roots,i), i);
		arr_print(stdout,&debug, "%.5lf\n" );
		putchar('\n');
//		convergence_study( roots.val[i], &debug );
	}
	

    plot_graph( f, x1, x2, &roots, starting_points );
	
//	print_table( &starting_points, &roots );

	Array *arr_to_free[]={&roots,&debug};
	arr_free_many(arr_to_free,2);
}

int main()
{
	Array starting_points = arr_linspace(1,7,7);
	
	test(f1, f1_der, 0.5, 10, &starting_points);	
	
	arr_free(&starting_points);
	
	return 0;
}

