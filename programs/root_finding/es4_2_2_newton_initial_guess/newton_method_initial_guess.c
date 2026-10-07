//Es15_2.c
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rootf.h"
#include "explot.h"
#include <string.h>

#define DOMAIN_PLOT_SIZE 1000
#define MANY_CHARS 2000
#define NTRIES 7

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

int advance_index(bool reset)
{
	static int i=0;
	if(reset)
		i=0;
	const int to_return = i;
	i++;
	return to_return;
}

void plot_graph( ParamFuncPtr f, ParamFuncPtr der,
	double x1, double x2, const Array *roots,
	const Array *initial_guesses, bool show_derivative[], int nroots,
		bool bracket)
{	
	SeriesSpec *series;
	SAFE_ALLOC(series, 5 + nroots);
	
	double x_line[] = {x1, x2};
	double x_line_zeroes[2] = {0};
	series[advance_index(true)] = eplot_init_series(x_line, x_line_zeroes, 2, 
		"lines", NULL, "rgb \"#50000000\" lw 1.2");
	
	
	Array *domain_derivatives;
	SAFE_ALLOC(domain_derivatives, nroots);
	Array *y_value_derivatives;
	SAFE_ALLOC(y_value_derivatives, nroots);
	Array x_init_guesses={0};
	Array y_init_guesses={0};
	
	int count_derivatives_to_show = 0;
	
	for(int i=0; i<nroots; i++)
		if(show_derivative)
		if(show_derivative[i])
			{
				Array *curr_domain =  domain_derivatives+i;
				Array *curr_yvalues = y_value_derivatives+i;
				
				const int pt_to_draw = 10;
				
				arr_null(curr_yvalues);
				arr_init(curr_yvalues, pt_to_draw);
				
				const double iguess = ARRP(initial_guesses, i);
				const double figuess = f(iguess, NULL);
				const double derivative = der(iguess, NULL);				
				const double xintersection = iguess -figuess / derivative;
				
				arr_append(&x_init_guesses, iguess);
				arr_append(&y_init_guesses, figuess);
				
				printf("%lf %lf %lf %lf\n", iguess, figuess, derivative,
					xintersection);
				
				*curr_domain = arr_linspace( fmin(iguess, xintersection)-0.3, 
					fmax(iguess, xintersection)+0.3, pt_to_draw );
				
				
				for(int j=0; j<pt_to_draw; j++)
					printf( "%lf\n",ARRP(curr_yvalues, j) = derivative*
						(ARRP(curr_domain,j) - iguess) + figuess);
				
				series[advance_index(false)] = 
					eplot_init_series(curr_domain->data, 
					curr_yvalues->data, pt_to_draw, "l dt 2", NULL, NULL);
				count_derivatives_to_show++;
			}
	
	
	Array domain = arr_linspace(x1,x2,DOMAIN_PLOT_SIZE);
	Array y_values = arr_map_par(&domain,f,NULL);

	series[advance_index(false)] = eplot_init_series(
		domain.data,y_values.data,domain.size, "l","f(x)",NULL);
	
	Array y_roots = arr_map_par(roots,f,NULL);
	
	series[advance_index(false)] = eplot_init_series(
		roots->data,y_roots.data,roots->size, "p pt 2","Roots",NULL);
	
	Array initial_guesses_y = arr_new(initial_guesses->size);
	
	if(!bracket)
	{
		arr_setv(&initial_guesses_y, -1);
		
		series[advance_index(false)] = eplot_init_series(
			initial_guesses->data, initial_guesses_y.data, initial_guesses->size,
			"p pt 4","Initial guesses", NULL);
	}else{
		arr_setv(&initial_guesses_y, 0);
		
		series[advance_index(false)] = eplot_init_series(
			initial_guesses->data, initial_guesses_y.data, initial_guesses->size,
			"p pt 5", "Interval delimiters", NULL);
	}
	
	
	char add_set[MANY_CHARS];
	int str_end = 0;
	
	str_end+=snprintf(add_set, sizeof add_set - str_end, 
	"set xrange [%.2lf: %.2lf]\n", x1, x2);
	if(!bracket)
		for(int i=0; i<initial_guesses->size; i++)
			str_end+=sprintf(add_set+str_end,"set arrow from %lf,%lf to %lf,"
			"%lf head filled lw 2.5 lc rgb \"#D45D5D\" front\n",
			ARRP(initial_guesses,i),ARR(initial_guesses_y,i),
			ARRP(roots,i),0.0);
	
	if(show_derivative)
		series[advance_index(false)] = eplot_init_series( x_init_guesses.data,
			y_init_guesses.data, count_derivatives_to_show, "p", "", NULL);
	
	Global gb_settings = eplot_init_global("Plot of function x^{-2} - sin(x)",
						"x","y",NULL,"newt_root_discussion","newt_root_discussion",NULL,
						add_set);
	
	
	
	eplot_multi_options(series, count_derivatives_to_show+(show_derivative ? 5 : 4),
		&gb_settings, true);
	
	
	
	arr_free(&domain);
	arr_free(&y_values);
	arr_free(&y_roots);
	arr_free(&initial_guesses_y);
	
	for (int i = 0; i < nroots; i++) {
		arr_free(domain_derivatives+i);
		arr_free(y_value_derivatives+i);
	}
	SAFE_FREE(domain_derivatives);
	SAFE_FREE(y_value_derivatives);
	SAFE_FREE(series);
	arr_free(&x_init_guesses);
	arr_free(&y_init_guesses);
	
}

void test( ParamFuncPtr f, ParamFuncPtr der, 
			double x1, double x2, const Array *starting_points,
			bool use_bracketing,bool show_derivative[], int nroots_searched)
{
	Array debug = {0};

	UNUSED(x1);
	UNUSED(x2);
	
	Array roots = arr_new(nroots_searched);
	
	for( int i = 0; i<nroots_searched; i++)
	{
		if(use_bracketing)
		{
			ARR(roots,i) = root_newt_bracket( f, NULL, der, NULL, ARRP(starting_points,2*i),
						ARRP(starting_points, 2*i+1),
						DEF_TOL, DEF_TOL, &debug );
		}
		else{
			ARR(roots,i) = root_newt( f, NULL, der, NULL, ARRP(starting_points,i),
						DEF_TOL, DEF_TOL, &debug );
		}

		printf("Intermediate values for finding the root %lf of index %d:\n", 
				ARR(roots,i), i);
		arr_print(stdout, &debug, "%.20lf\n" );
		putchar('\n');
	}

    plot_graph( f, der, x1, x2, &roots, starting_points, show_derivative, 
		nroots_searched, use_bracketing );
	
//	print_table( &starting_points, &roots );

	Array *arr_to_free[]={&roots,&debug};
	arr_free_many(arr_to_free,2);
}

int main()
{
	Array starting_points = arr_linspace(1, NTRIES, NTRIES);
	bool show_derivative[NTRIES]={false};
	show_derivative[1]=true;
	show_derivative[4]=true;
	
	test(f1, f1_der, 0.5, 10, &starting_points, false, show_derivative, NTRIES);
	
	double starting_points_bracketing_carr[] = {0.5, 2, 2, 5, 5, 8, 8, 10};
	Array starting_points_bracketing = 
		arr_asarr(starting_points_bracketing_carr, 8);
	
	test(f1, f1_der, 0.5, 10, &starting_points_bracketing, true, NULL, 4);
	
	arr_free(&starting_points);
	
	return 0;
}

