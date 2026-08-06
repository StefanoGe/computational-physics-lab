#include "interpolation.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "explot.h"

#define PLOT_DOMAIN_LENGTH 100

MAKE_PARAM_FUNC_WRAP(log)
MAKE_PARAM_FUNC_WRAP(fabs)
MAKE_PARAM_FUNC_WRAP(cosh)
MAKE_PARAM_FUNC_WRAP(tanh)

Array barf_fit(const BarFit *barf, SeriesSpec *series,
			int *added_series, const Array *domain)
{
	Array approx_pol = arr_new(PLOT_DOMAIN_LENGTH);
	for(int i=0; i<PLOT_DOMAIN_LENGTH; i++)
		ARR(approx_pol, i)=interp_barf_get_value(barf,ARRP(domain,i));
		
	arr_print_inline(&approx_pol,"%lf ",true,true);
	
	series[*added_series] = eplot_init_series(domain->data, 
	approx_pol.data, PLOT_DOMAIN_LENGTH, "l", "approx pol", NULL);
	
	(*added_series)++;
	
	return approx_pol;
}

void test_int_eq_spaced(ParamFuncPtr f, double x1, double x2, int nnodes[],
						int nruns)
{	
	ParamFunc fp = {f,NULL};
	BarFit barf = {0};
	interp_barf_add_func(&barf,&fp);
	
	SeriesSpec *series;
	SAFE_ALLOC(series,nruns+1);
	int added_series = 1;
	
	Array *approx_pols;
	SAFE_ALLOC(approx_pols,nruns);
	Array domain = arr_linspace(x1,x2,PLOT_DOMAIN_LENGTH);
	Array true_values = arr_map_par(&domain,f,NULL);
	arr_print_inline(&true_values,"%lf ",true,true);
	
	Global gb_settings = eplot_init_global(
		"Interpolation with equally spaced nodes", "x", "y", NULL, NULL, 
		"equally_spaced_barycentric_fit", "top left", NULL);
		
	series[0] = eplot_init_series(
	domain.data, true_values.data, PLOT_DOMAIN_LENGTH, "l", "true values", NULL);
	
	for(int i=0; i<nruns; i++)
	{
		interp_barf_init_eq_points(&barf,x1, x2, nnodes[i]);
		interp_barf_add_fvalues(&barf,NULL);
		approx_pols[i]=barf_fit(&barf,series, &added_series,&domain);
	}
	
	eplot_multi(series, nruns+1, &gb_settings);
	
	interp_barf_free(&barf);
	arr_free(&domain);
	for(int i=0;i<nruns;i++)
		arr_free(approx_pols+i);
	SAFE_FREE(series);
	SAFE_FREE(approx_pols);
}

int main()
{
	int n_nodes_1 []= {2, 3, 4};
	test_int_eq_spaced( log_wrap, 1,10,n_nodes_1, 3 );
	test_int_eq_spaced( tanh_wrap, -3, 2, n_nodes_1, 3 );
	test_int_eq_spaced( cosh_wrap, -1, 3, n_nodes_1,3 );
	int n_nodes_2[] = {3, 5, 7};
	test_int_eq_spaced( fabs_wrap, -2, 1, n_nodes_2,3 );
	
	exit(EXIT_SUCCESS);
}
