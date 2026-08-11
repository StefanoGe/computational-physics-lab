#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "interpolation.h"
#include "explot.h"

#define N_DOMAIN 1000

double f1(double x, void *null)
{
	(void)null;
	return cosh(sin(x));
}

void test_cheb(ParamFuncPtr f, double x1, double x2, int n_nodes)
{
	ParamFunc fp = {f,NULL};
	BarFit barf = interp_barf_new_cheb2_points(x1,x2,n_nodes);
	interp_barf_add_fvalues(&barf,&fp);
	
	Array domain = arr_linspace(x1,x2,N_DOMAIN);
	Array approx_pol = arr_map_par(&domain,interp_barf_get_value_wrap,&barf);
	arr_print_inline(&approx_pol,"%lf ",true,true);
	Array true_values = arr_map_par(&domain, f, NULL);
	
	Global gb_settings = eplot_init_global("Polynomial interpolation of "
				"cosh(sin x)) with Chebyshev nodes","x","y",NULL,
				"cheb_example",
				NULL,NULL,NULL);
	
	SeriesSpec data = eplot_init_series(domain.data,true_values.data,
			N_DOMAIN,"l","cosh(sin(x))",NULL);
	
	SeriesSpec plot_fit = eplot_init_series(domain.data,approx_pol.data,
			N_DOMAIN,"l dt 3","Interpolant polynomial",NULL);
			
//	Array node_indices = arr_linspace(0,n_nodes-1,1);
			
	SeriesSpec nodes = eplot_init_series(barf.points,barf.f_values,
		n_nodes,"p","Nodes",NULL);
			
	SeriesSpec series[]={data,plot_fit,nodes};
	
	eplot_multi(series,3,&gb_settings);
	
	Array *arr_to_free[] = {&domain,&approx_pol,&true_values};
	arr_free_many(arr_to_free,3);
	interp_barf_free(&barf);
}

int main()
{
	test_cheb(f1,0,2*PI,40);
	
	exit(EXIT_SUCCESS);
}

