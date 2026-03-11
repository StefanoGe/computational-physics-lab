#include <stdlib.h>
#include <stdio.h>
#include "linearsys.h"
#include <math.h>

MAKE_PARAM_FUNC_WRAP(log)

double const_1( double x, void *p )
{
	UNUSED(p);
	UNUSED(x);
	return 1;
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
	
	ParamFunc f_log=param_func_null(log_wrap);
	ParamFunc f_const_1=param_func_null(const_1);
	
	Array log_periods = arr_new( 8 );
	for( int i = 0; i < 8; i++)
		ARR(log_periods,i) = log(periods[i]);
	
	LinearModel model = linear_model_alloc( 2 );
	model.funcs[0] = f_const_1;
	model.funcs[1] = f_log;
	
	Array coeffs={0};
	linst_lsqr_fit_linear( &distances_arr, &log_periods, &model,
		&coeffs);
	
	arr_print_inline( &coeffs, "%g", true, true );
	
	linear_model_free(&model);
	arr_free(&distances_arr);
	arr_free(&coeffs);
	arr_free(&log_periods);
	
}

int main()
{
	kepler_law_int();
	
	return 0;
}
