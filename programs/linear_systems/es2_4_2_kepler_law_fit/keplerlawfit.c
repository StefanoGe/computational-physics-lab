#include <stdlib.h>
#include <stdio.h>
#include "linearsys.h"
#include <math.h>
#include "explot.h"

MAKE_PARAM_FUNC_WRAP(log)

double const_1( double x, void *p )
{
	UNUSED(p);
	UNUSED(x);
	return 1;
}

void kepler_plot( double distances[], double periods[],
	double alpha, double c)
{
	Global gb_settinga = eplot_init_global("Kepler law", "Distances from the sun [Mkm]",
		"Orbital periods [days]", NULL, NULL, "kepler_plot", NULL, NULL);
	
	SeriesSpec exp_data = eplot_init_series(distances, periods, 8, 
		"p", "Data", NULL);
		
	int count_data = 1000;
		
	Array distances_calc = arr_linspace(0, 5000, count_data);
	Array periods_calc = arr_new(count_data);
	
	for(int i=0; i<count_data; i++)
		ARR(periods_calc,i) = c * pow( ARR(distances_calc,i), alpha );
		
	SeriesSpec calc_data = eplot_init_series(distances_calc.data, 
		periods_calc.data, count_data, "l", "Fit", NULL);
	
	SeriesSpec series[] = {exp_data, calc_data};
	
	eplot_multi(series, 2, &gb_settinga);
	
	arr_free(&distances_calc);
	arr_free(&periods_calc);
}

void kepler_law_int()
{
	double distances[] = {
		57.59,
		108.11,
		149.57,
		227.84,
		778.14,
		1427.0,
		2870.3,
		4499.9};
		
	Array distances_arr = arr_asarr(distances, 8);
		
	double periods []= {
		87.99,
		224.7,
		365.26,
		686.98,
		4332.4,
		10759.0,
		30684.0,
		60188.0};
		
	Array periods_arr = arr_asarr(periods, 8);
	
	ParamFunc f_log=param_func_null(log_wrap);
	ParamFunc f_const_1=param_func_null(const_1);
	
	Array log_periods = arr_map(&periods_arr, log);

	LinearModel model = linear_model_alloc( 2 );
	model.funcs[0] = f_const_1;
	model.funcs[1] = f_log;
	
	Array coeffs={0};
	linst_lsqr_fit_linear_lup( &distances_arr, &log_periods, &model, &coeffs);
	
	const double alpha = ARR(coeffs,1);
	const double c = exp(ARR(coeffs,0));

	
	printf("Found:\nalpha = %g\nc = %g\n\n", alpha, c);
	fflush(stdout);
	kepler_plot(distances, periods, alpha, c);
	
	linear_model_free(&model);
	arr_free(&distances_arr);
	arr_free(&coeffs);
	arr_free(&log_periods);
	arr_free(&periods_arr);
}

int main()
{
	kepler_law_int();
	
	return 0;
}
