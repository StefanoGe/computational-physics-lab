#include "interpolation.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "explot.h"

#define PLOT_DOMAIN_LENGTH 100
#define MAX_NRUNS 10
#define MAX_POINTTYPE_CHARS 20
#define MAX_NAMEPOL_CHARS 50

MAKE_PARAM_FUNC_WRAP(log)
MAKE_PARAM_FUNC_WRAP(fabs)
MAKE_PARAM_FUNC_WRAP(cosh)
MAKE_PARAM_FUNC_WRAP(tanh)

Array barf_fit(const BarFit *barf, SeriesSpec *series,
			int *added_series, const Array *domain, char *pointtype, int index,
			char *polynomial_plot_name)
{
	Array approx_pol = arr_new(PLOT_DOMAIN_LENGTH);
	for(int i=0; i<PLOT_DOMAIN_LENGTH; i++)
		ARR(approx_pol, i)=interp_barf_get_value(barf,ARRP(domain,i));
		
//	arr_print_inline(&approx_pol,"%lf ",true,true);
	
	sprintf(polynomial_plot_name, "%d nodes", barf->npoints);
	
	series[*added_series] = eplot_init_series(domain->data, 
	approx_pol.data, PLOT_DOMAIN_LENGTH, "l", polynomial_plot_name, NULL);
	
	int pt = (index+1>2 ? index+2 : index+1);
	
	sprintf(pointtype,"p pt %d", pt);
	
	series[*added_series+1] = eplot_init_series(barf->points, barf->f_values,
		barf->npoints, pointtype, NULL, NULL);
	
	(*added_series)+=2;
	
	return approx_pol;
}

void test_int_eq_spaced(ParamFuncPtr f, double x1, double x2, int nnodes[],
						int nruns, const char *fname)
{	
	ParamFunc fp = {f,NULL};
	BarFit *barfs;
	SAFE_ALLOC(barfs,nruns);
	for(int i=0; i<nruns; i++)
		interp_barf_add_func(barfs+i,&fp);
	
	SeriesSpec *series;
	SAFE_ALLOC(series,2*nruns+1);
	int added_series = 1;
	
	Array *approx_pols;
	SAFE_ALLOC(approx_pols,nruns);
	Array domain = arr_linspace(x1,x2,PLOT_DOMAIN_LENGTH);
	Array true_values = arr_map_par(&domain,f,NULL);
//	arr_print_inline(&true_values,"%lf ",true,true);
	
	char pointtype[MAX_NRUNS][MAX_POINTTYPE_CHARS];
	char polynomial_plot_name[MAX_NRUNS][MAX_NAMEPOL_CHARS];
	char title[100];
	sprintf(title,"Interpolation with equally spaced nodes of %s",fname);
	
	Global gb_settings = eplot_init_global( 
	title,"x", "y", NULL, NULL, 
		"equally_spaced_barycentric_fit", "top left", NULL);
		
	series[0] = eplot_init_series(
	domain.data, true_values.data, PLOT_DOMAIN_LENGTH, "l", "Exact function", NULL);
	
	for(int i=0; i<nruns; i++)
	{
		interp_barf_init_eq_points(barfs+i,x1, x2, nnodes[i]);
		interp_barf_add_fvalues(barfs+i,NULL);
		approx_pols[i]=barf_fit(barfs+i,series, &added_series,&domain,
			pointtype[i], i,polynomial_plot_name[i]);
	}
	
	eplot_multi(series, 2*nruns+1, &gb_settings);
	
	for(int i=0; i<nruns; i++)
		interp_barf_free(barfs+i);
	SAFE_FREE(barfs);
	arr_free(&domain);
	arr_free(&true_values);
	for(int i=0;i<nruns;i++)
		arr_free(approx_pols+i);
	SAFE_FREE(series);
	SAFE_FREE(approx_pols);
}

int main()
{
	int n_nodes_1 []= {2, 3, 4};
	test_int_eq_spaced( log_wrap, 1, 10,n_nodes_1, 3, "log" );
	test_int_eq_spaced( tanh_wrap, -3, 2, n_nodes_1, 3, "tanh" );
	test_int_eq_spaced( cosh_wrap, -1, 3, n_nodes_1, 3, "cosh" );
	int n_nodes_2[] = {3, 5, 7};
	test_int_eq_spaced( fabs_wrap, -2, 1, n_nodes_2,3, "|.|" );
	
	exit(EXIT_SUCCESS);
}
