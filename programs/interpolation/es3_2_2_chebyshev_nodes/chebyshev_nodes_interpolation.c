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

	double maximum_difference(double *arr1, double *arr2, int size)
	{
		double maximum = 0;
		for(int i=0; i<size; i++)
			maximum = fmax(maximum, fabs(arr1[i]-arr2[i]));
		return maximum;	
	}

void test_cheb(ParamFuncPtr f, double x1, double x2, int n_nodes)
{
	ParamFunc fp = {f,NULL};
	BarFit barf = interp_barf_new_cheb2_points(x1,x2,n_nodes);
	interp_barf_add_fvalues(&barf,&fp);
	
	double *generic_formula_weights = interp_build_weights(barf.points, 
		barf.npoints);

	// Since weights have a different constant factor preceding them
	// in the generic and analytic formula, we must scale in order to
	// compare them

	const double scale = generic_formula_weights[2];
	
	for(int i=0; i<n_nodes; i++)
	{
		generic_formula_weights[i] /= scale;
//		fprintf(stderr, "%g %g\n", barf.weights[i], generic_formula_weights[i]);
	}
	
	fprintf(stderr, ANSI_GREEN"In the evaluation of the weight for function "
		ANSI_RED"cosh(sin(x)) "ANSI_GREEN"the maximum \ndifference between "
		"the weights of the general formlula and the specialized one \n"
		"for equally spaced nodes is: "ANSI_RESET"%g\n", 
			maximum_difference(barf.weights, generic_formula_weights, n_nodes));
	
	Array domain = arr_linspace(x1,x2,N_DOMAIN);
	Array approx_pol = arr_map_par(&domain,interp_barf_get_value_wrap,&barf);
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
	
	eplot_multi_options(series,3,&gb_settings,false);
	
	Array *arr_to_free[] = {&domain,&approx_pol,&true_values};
	arr_free_many(arr_to_free,3);
	interp_barf_free(&barf);
	SAFE_FREE(generic_formula_weights);
}

int main()
{
	test_cheb(f1,0,2*PI,40);
	
	exit(EXIT_SUCCESS);
}

