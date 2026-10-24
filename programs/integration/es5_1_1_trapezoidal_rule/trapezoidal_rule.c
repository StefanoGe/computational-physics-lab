// Test of trapezoidal rule of integration and study of its convergence
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "integration.h"
#include "explot.h"
#include "genutil.h"
#include "linearsys.h"

#define MAX_NODE_POWER 10
#define NTESTS 6

double id(double x,void *unused)
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

double power_minus_two(double x)
{
	return 1/x/x;
}

double f1( double x, void * unused )
{
	UNUSED(unused);
	return x * log(1 + x);
}

double f2( double x, void * unused )
{
	UNUSED(unused);
	return x*x*atan(x);
}

double f3( double x, void * unused )
{
	UNUSED(unused);
	return exp(x) * cos(x);
}

double f4( double x, void * unused )
{
	UNUSED(unused);
	const double A = sqrt( 2 + x * x );
	return atan( A ) / ( 1 + x * x )/A;
}

double f5( double x, void * unused )
{
	UNUSED(unused);
	if (x == 0)
		return 0;
	return sqrt(x) * log(x);
}

double f6( double x, void * unused )
{
	UNUSED(unused);
	return sqrt( 1 - x * x );
}

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


void plot_graph(const Array *errors,char *fnames[])
{
	Array n_node_list = arr_new(MAX_NODE_POWER);
	for( int i = 1; i <= MAX_NODE_POWER; i++ )
		ARR(n_node_list,i-1) = 10 * pown(2, i);
	Array y_second_order_reference_line = arr_map(&n_node_list, power_minus_two);
	
	Global gb_settings = eplot_init_global("Trapezoidal rule error plot",
		"Number of nodes","Error","xy","trap_err","trap_err","bottom left",NULL);
	
	SeriesSpec *series;
	SAFE_ALLOC(series, NTESTS+1);
	
	for(int i=0; i<NTESTS; i++)
		series[i] = eplot_init_series(n_node_list.data,
			(errors+i)->data,(errors+i)->size,"lp",fnames[i],NULL);
	
	series[NTESTS] = eplot_init_series(n_node_list.data,
		y_second_order_reference_line.data,n_node_list.size,
			"l" ,"1/x^2 ref line",NULL);

	eplot_multi_options(series, NTESTS+1, &gb_settings, false);
	
	SAFE_FREE(series);
	arr_free(&n_node_list);
}

void integration_study( const ParamFuncPtr f, double x1, double x2, 
			double true_value, Array *errors, int f_index)
{
	printf(
			ANSI_GREEN "Study of convergence of the integral of"
			ANSI_RED " function %d \n"
			ANSI_YELLOW "between endpoints %.2lf and %.2lf\n"
			"True value: %g\n\n" ANSI_RESET, f_index,x1,x2,true_value);
	printf("n_nodes | Error\n");
	for( int i = 1; i <= MAX_NODE_POWER; i++ )
	{
		int n_nodes = 10 * pow_int(2, i);
		double integral = intg_trap( f, NULL, x1, x2, n_nodes );
		double error = true_value - integral;
		//Absolute error or no? Yes for the log plot
		arr_append( errors, fabs(error) );
		printf( "%6d  | %.2e\n", n_nodes, error );
	}

	Array log_nodes = arr_linspace(1, MAX_NODE_POWER, MAX_NODE_POWER);
	Array log_errs = arr_map(errors,log2);
	Array coeffs_fit = {0};
	
	ParamFuncPtr fs[] = {id,one};
	linst_lsqr_fit_linear_qr_wrap(&log_nodes,&log_errs,fs,NULL,2,&coeffs_fit);
	
	printf( ANSI_YELLOW"\nParameters of fitted line in log-scale:\n"
			ANSI_RESET"m = %.3lf | q = %.3lf\n",
			ARR(coeffs_fit,0),
			ARR(coeffs_fit,1));
	
	puts( "\n----------------------------\n");
	
	arr_free(&log_nodes);
	arr_free(&log_errs);
	arr_free(&coeffs_fit);
}

int main()
{	
	Array *errors;
	SAFE_ALLOC(errors,NTESTS);
	for(int i=0; i<NTESTS; i++)
		arr_null(errors+i);

	puts("----------------------------\n");	
	integration_study( f1, 0, 1, 0.25,errors,1);
	integration_study( f2, 0, 1, (PI - 2 + 2 * log(2)) / 12, errors+1,2);	
	integration_study( f3, 0, PI/2, ( exp( PI/2 ) - 1 ) / 2,errors+2,3);
	integration_study( f4, 0, 1, 5 * PI*PI/96,errors+3,4);
	integration_study( f5, 0, 1, -4.0/9,errors+4,5);
	integration_study( f6, 0, 1, PI/4,errors+5,6);
	
	
	char *fnames[]={"f_1(x)","f_2(x)","f_3(x)",
			"f_4(x)","f_5(x)","f_6(x)"};
	
	plot_graph(errors,fnames);
	
	for(int i=0; i<NTESTS; i++)
		arr_free(errors+i);
	
	exit(EXIT_SUCCESS);
}

