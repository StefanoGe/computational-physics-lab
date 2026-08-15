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
	printf("%lf\n", root_bis( f1, NULL, -2.5, 0, 1e-10 , NULL) );
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

void test2(double x1, double x2)
{
	Array debug = {0};
	arr_reserve(&debug,1);
	
	const double root = root_bis( f2, NULL, x1, x2, 0, &debug );
	printf( "Found root in x = %.20lf\n", root );
	
//	arr_print(stdout,&debug,"%.20g\n");
	
	Array log_abs_err = arr_new(debug.size-1);
	for(int i=0; i<debug.size-1; i++)
		ARR(log_abs_err,i) = -log10 ( fabs( ARR(debug,i) - root ) );
		
	Array indices = arr_linspace(1,debug.size-1,debug.size);
	
	Array fit_result = arr_new(2);
	
	ParamFuncPtr funcs_for_fit[]={id,one};
	
//	arr_print_inline(&log_abs_err,"%lf ",true,true);
	
	linst_lsqr_fit_linear_qr_wrap(&indices,&log_abs_err,funcs_for_fit,NULL,
				2,&fit_result);
	
//	printf("%lf %lf\n", pow(10,ARR(fit_result,0)), ARR(fit_result,1));
	
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
	
	eplot_multi(series,3,&gb_setting);
	
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

// y = - log10(|diff|)

// se relazione è |diff| = (1/2)^n

// allora ho y = - log10((1/2)^n) = - n * log10(1/2) = n * log10( 2 )

