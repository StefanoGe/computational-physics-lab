#include "linearsys.h"
#include "genutil.h"
#include "explot.h"
#include <stdio.h>
#include <math.h>

#define N_XS 100
#define NORMALIZATION_CONSTANT 2006.787453080206
#define LUP 0
#define MGS 1
#define QLESS 2

double monomyal(double x, void *exp){return pow(x, *(int*)exp);}

double h( double t )
{
	return exp( sin( 4 * t ) );
}

double h_norm(double t)
{
	return h(t)/NORMALIZATION_CONSTANT;
}

void plot_fit( const Array *coeffs, const Array *x_data,
	const Array *y_data, const LinearModel *model, int n_data)
{
	
	const int npoints = 500;
	Array x_data_plot = arr_linspace(0, 1, npoints);
	Array y_data_plot = arr_new(npoints);
	
	for(int i=0; i<npoints; i++)
		ARR(y_data_plot, i) = linst_linear_eval(ARR(x_data_plot,i), model,
			coeffs);

	Global gb_settings = eplot_init_global(
		"Interpolation of exp(sin(4t)) with monomyals", 
		"t",
		"y", NULL, NULL, "expsin_pol_qr", NULL, NULL);
	
	SeriesSpec fit_data = eplot_init_series(x_data_plot.data, y_data_plot.data, npoints, 
		"l", "Fit", NULL);
		
		
	SeriesSpec true_data = eplot_init_series(x_data->data, y_data->data, n_data, 
		"p", "Exact function", NULL);
		
	SeriesSpec series[] = {fit_data, true_data};
	
	eplot_multi(series, 2, &gb_settings);
	
	
	arr_free(&x_data_plot);
	arr_free(&y_data_plot);
}

void fit_exp_sin(int n_pol, int mode, bool do_plot)
{
	const int n_data = 100;
	Array x_data = arr_linspace( 0, 1, n_data );
	Array y_data = arr_map( &x_data, h_norm );
	Array coeffs = {0};
	
	int *exponents = malloc(n_pol*sizeof(int));
	
	LinearModel model = linear_model_alloc(n_pol);
	for(int i=0; i<n_pol; i++)
	{
		exponents[i]=i;
		model.funcs[i].func=monomyal;
		model.funcs[i].params=exponents+i;
	}
	
	if(mode==LUP)
		linst_lsqr_fit_linear_lup(&x_data, &y_data, &model, &coeffs);
	else if(mode==MGS)
		linst_lsqr_fit_linear_qr(&x_data, &y_data, &model, &coeffs);
	else if(mode==QLESS)
		linst_lsqr_fit_linear_qless(&x_data, &y_data, &model, &coeffs);
	else
		raiseErr("Select an appropriate mode");

	char *used_decomposition = (mode==0 ? "LUP" : (mode==1 ? "QR" : "Q-less") );
	printf("Found coeffs using %s decomposition:\n", used_decomposition);
	arr_print_inline(&coeffs, "%g", true, true);

	printf("Last coefficient c%d: %.10g\n\n", n_pol, ARR(coeffs, n_pol-1));

	if(do_plot)
		plot_fit( &coeffs, &x_data, &y_data, &model, n_data);

	Array *arr_to_free[]={&x_data, &y_data, &coeffs};
	arr_free_many(arr_to_free, 3);
	free(exponents);
	linear_model_free(&model);
}

int main()
{
	fit_exp_sin(15, LUP, false);
	fit_exp_sin(15, MGS, false);
	fit_exp_sin(15, QLESS, true);
	
	exit(EXIT_SUCCESS);
}
