#include <stdio.h>
#include "rootf.h"
#include <math.h>
#include "explot.h"
#include "linearsys.h"

#define DOMAIN_PLOT_SIZE 1000

double f1( double x, void * unused )
{
	UNUSED(unused);
	return x + 1;
}

void test1()
{	
	printf(ANSI_GREEN"Testing the bisection method rootfinder on the function \n"
		"f(x) = "ANSI_RED"x+1"ANSI_GREEN
		", we find a root in x = "ANSI_YELLOW"%lf\n"ANSI_RESET"\n"
		"--------------------------------------------------------\n\n", 
		root_bis( f1, NULL, -2.5, 0, 1e-10 , NULL) );
}

double f2(double x, void *unused)
{
	UNUSED(unused);
	return x*x*x - 7*x*x + 14*x -6;
}

double id(double x, void *unused)
{
	UNUSED(unused);
	return x;
}

double one(double x, void *unused)
{
	UNUSED(unused);
	UNUSED(x);
	return 1;
}

double line(double x, void *params)
{
	return x*((double *)params)[0] + ((double *)params)[1];
}

void plot_func(double root)
{
	Array domain = arr_linspace(0,1,DOMAIN_PLOT_SIZE);
	Array y_values = arr_map_par(&domain, f2, NULL);
	double zero = 0;
	
	Global gb_settings = eplot_init_global("Plot of x^3-7x^2+14x-6",
		"x", "y", NULL, "bisection_sample_plot", "", NULL,NULL);
		
	SeriesSpec fnc = eplot_init_series(domain.data, y_values.data,
		DOMAIN_PLOT_SIZE,"l","f(x)",NULL);
		
	SeriesSpec root_point = eplot_init_series(&root, &zero,
		1,"p","f(x)",NULL);
		
	SeriesSpec series[] = {fnc, root_point};
		
	eplot_multi_options(series, 2, &gb_settings, false);
	
	arr_free(&y_values);
	arr_free(&domain);
}

void test2(double x1, double x2)
{
	Array debug = {0};
	arr_reserve(&debug,1);
	
	const double root = root_bis( f2, NULL, x1, x2, 0, &debug );
	printf(ANSI_GREEN"Testing the bisection method rootfinder on the function \n"
		"f(x) = "ANSI_RED"x^3-7x^2+14x-6"ANSI_GREEN 
		", considering only the interval [0,1]\nwe find a root in x = "
		ANSI_YELLOW"%lf\n"ANSI_RESET, 
		root );
	
//	arr_print(stdout,&debug,"%.20g\n");
	
	Array log_abs_err = arr_new(debug.size-1);
	for(int i=0; i<debug.size-1; i++)
		ARR(log_abs_err,i) = -log10(fabs(ARR(debug,i)-root));
		
	Array indices = arr_linspace(1,debug.size-1,debug.size-1);
	
	Array fit_result = arr_new(2);
	
	ParamFuncPtr funcs_for_fit[]={id,one};
	
//	arr_print_inline(&log_abs_err,"%lf ",true,true);
	
	linst_lsqr_fit_linear_qr_wrap(&indices, &log_abs_err,funcs_for_fit,NULL,
				2,&fit_result);
	
	printf(
		ANSI_GREEN"Slope of the line is: "ANSI_YELLOW"%g\n"
		ANSI_GREEN"The estimated asymptotic error constant is: "
		ANSI_YELLOW"%g\n"ANSI_RESET"\n"
		"--------------------------------------------------------\n\n", 
		pow(10,ARR(fit_result,0)),
		pow(10,-ARR(fit_result,0)));
	
	plot_func(root);
	
	Array fitted_line = arr_map_par(&indices,line,fit_result.data);
	
	double line_params[]={log10(2),0};
	
	Array exp_line_data = arr_map_par(&indices,line,line_params);
	
	Global gb_setting = eplot_init_global("Convergence of bisection method",
		"n","-log_{10}|x_n - x_{est}|",NULL,"convergence_bis",NULL,"bottom right",NULL);
	
	SeriesSpec conv = eplot_init_series(indices.data,log_abs_err.data,
		debug.size-1,"lp","Error found at each step",NULL);
		
	SeriesSpec fit_line = eplot_init_series(indices.data, fitted_line.data,
		debug.size-1,"l","Line fit with all points",NULL);
		
	SeriesSpec expected_line = eplot_init_series(indices.data, exp_line_data.data,
		debug.size-1,"l","Upper bound error",NULL);
		
	SeriesSpec series[]={conv,fit_line,expected_line};
	
	putchar('\n');
	
	eplot_multi_options(series, 3, &gb_setting, false);
	
	arr_free(&debug);
	arr_free(&log_abs_err);
	arr_free(&indices);
	arr_free(&fit_result);
}



int main()
{
	test1();
	test2(0,1);
/*
	test2(0.4, 1);
	test2(0.25, 0.66);
*/
	
	exit(EXIT_SUCCESS);
}
