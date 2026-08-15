//Es15_2.c
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rootf.h"
#include "explot.h"
#include <string.h>

#define DOMAIN_PLOT_SIZE 1000
#define MANY_CHARS 2000

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

void plot_graph( ParamFuncPtr f, double x1, double x2, const Array *roots,
	const Array *initial_guesses, bool show_derivative[])
{
	const int nroots = roots->size;
	
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
	
	SeriesSpec series[]={graph,root_rep,init_guesses_series};
	
	char set_arrows[MANY_CHARS];
	int str_end = 0;
	
	for(int i=0; i<initial_guesses->size; i++)
		str_end+=sprintf(set_arrows+str_end,"set arrow from %lf,%lf to %lf,"
		"%lf head filled lw 2.5 lc rgb \"#D45D5D\"\n",ARRP(initial_guesses,i),ARR(initial_guesses_y,i),
		ARRP(roots,i),0.0);
	
	Global gb_settings = eplot_init_global("Plot of function x^{-2} - sin(x)",
						"x","y",NULL,"newt_root_discussion",NULL,NULL,
						set_arrows);
	
	int count_derivatives_to_show = 0;
	for(int i=0; i<nroots; i++)
		if(show_derivative[i])
			count_derivatives_to_show++;
			
	
	
	eplot_multi(series,3,&gb_settings);
	
	arr_free(&domain);
	arr_free(&y_values);
	arr_free(&y_roots);
	arr_free(&initial_guesses_y);
}

void test( ParamFuncPtr f, ParamFuncPtr der, 
			double x1, double x2, const Array *starting_points,
			bool use_bracketing,bool show_derivative[])
{
	Array debug = {0};

	UNUSED(x1);
	UNUSED(x2);

	const int nroots_searched = starting_points->size;
	Array roots = arr_new(nroots_searched);
	
	for( int i = 0; i < nroots_searched; i++)
	{
		if(use_bracketing)
		{
			
		}
		else{
			ARR(roots,i) = root_newt( f, NULL,der,NULL,ARRP(starting_points,i),
						DEF_TOL, DEF_TOL, &debug );
		}
		printf("Intermediate values for finding the root %lf of index %d:\n", 
				ARR(roots,i), i);
		arr_print(stdout,&debug, "%.5lf\n" );
		putchar('\n');
//		convergence_study( roots.val[i], &debug );
	}
	

    plot_graph( f, x1, x2, &roots, starting_points, show_derivative );
	
//	print_table( &starting_points, &roots );

	Array *arr_to_free[]={&roots,&debug};
	arr_free_many(arr_to_free,2);
}

int main()
{
	Array starting_points = arr_linspace(1,7,7);
	bool show_derivative[]={0};
	show_derivative[1]=true;
	show_derivative[4]=true;
	
	test(f1, f1_der, 0.5, 10, &starting_points,false,show_derivative);
	
	arr_free(&starting_points);
	
	return 0;
}

