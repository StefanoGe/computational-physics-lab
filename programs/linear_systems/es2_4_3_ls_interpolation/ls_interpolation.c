#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "linearsys.h"
#include "explot.h"

double e_sin(double x)
{
	return exp(sin(x-1));
}

double const_1( double x, void *p )
{
	UNUSED(p);
	UNUSED(x);
	return 1;
}

double monomyal(double x, void *exp){return pow(x, *(int*)exp);}

double cospar(double x, void *k){return cos( *(int*)k * x);}

double sinpar(double x, void *k){return sin( *(int*)k * x);}

void approx_e_sin_pol(int n_pol)
{
	Array x_data = arr_linspace( 2 * PI / 60, 2*PI, 60 );
	Array y_data = arr_map( &x_data, e_sin );
	Array coeffs = {0};
	
	int *exponents = malloc(n_pol*sizeof(int));
	
	LinearModel model = linear_model_alloc(n_pol);
	for(int i=0; i<n_pol; i++)
	{
		exponents[i]=i;
		model.funcs[i].func=monomyal;
		model.funcs[i].params=exponents+i;
	}
	
	linst_lsqr_fit_linear_lup(&x_data, &y_data, &model, &coeffs);
	
	const int npoints = 200;
	Array x_data_plot = arr_linspace(0, 2*PI, npoints);
	Array y_data_plot = arr_new(npoints);
	
	for(int i=0; i<npoints; i++)
		ARR(y_data_plot, i) = linst_linear_eval(ARR(x_data_plot,i), &model,
			&coeffs);
	
	Global gb_settings = eplot_init_global(
		"Interpolation of exp(sin(t-1)) with monomyals", 
		"t",
		"y", NULL, NULL, "expsin_pol", NULL, NULL);
	
	SeriesSpec fit_data = eplot_init_series(x_data_plot.data, y_data_plot.data, npoints, 
		"l", "Fit", NULL);
		
		
	SeriesSpec true_data = eplot_init_series(x_data.data, y_data.data, 60, 
		"p", "Exact function", NULL);
		
	SeriesSpec series[] = {fit_data, true_data};
	
	eplot_multi(series, 2, &gb_settings);
	
	Array *arr_to_free[]={&x_data, &y_data, &coeffs, &x_data_plot, &y_data_plot};
	arr_free_many(arr_to_free, 5);
	free(exponents);
	linear_model_free(&model);
}

void approx_e_sin_periodic()
{
	Array x_data = arr_linspace( 2 * PI / 60, 2*PI, 60 );
	Array y_data = arr_map( &x_data, e_sin );
	Array coeffs = {0};
	
	LinearModel model = linear_model_alloc(5);
	int k[]={1,2};
	model.funcs[0].func=const_1;
	model.funcs[0].params=NULL;
	model.funcs[1].func=cospar;
	model.funcs[1].params=k;
	model.funcs[2].func=sinpar;
	model.funcs[2].params=k;
	model.funcs[3].func=cospar;
	model.funcs[3].params=k+1;
	model.funcs[4].func=sinpar;
	model.funcs[4].params=k+1;
	

	
	linst_lsqr_fit_linear_lup(&x_data, &y_data, &model, &coeffs);
	
	const int npoints = 200;
	Array x_data_plot = arr_linspace(0, 2*PI, npoints);
	Array y_data_plot = arr_new(npoints);
	
	for(int i=0; i<npoints; i++)
		ARR(y_data_plot, i) = linst_linear_eval(ARR(x_data_plot,i), &model,
			&coeffs);
	
	Global gb_settings = eplot_init_global(
		"Interpolation of exp(sin(t-1)) with periodic functions", 
		"t",
		"y", NULL, NULL, "expsin_per", NULL, NULL);
	
	SeriesSpec fit_data = eplot_init_series(x_data_plot.data, y_data_plot.data, npoints, 
		"l", "Fit", NULL);
		
		
	SeriesSpec true_data = eplot_init_series(x_data.data, y_data.data, 60, 
		"p", "Exact function", NULL);
		
	SeriesSpec series[] = {fit_data, true_data};
	
	eplot_multi(series, 2, &gb_settings);
	
	Array *arr_to_free[]={&x_data, &y_data, &coeffs, &x_data_plot, &y_data_plot};
	arr_free_many(arr_to_free, 5);
	linear_model_free(&model);
}

int main()
{

	approx_e_sin_pol(7);
	approx_e_sin_periodic();
	
	return 0;
}





