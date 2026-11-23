#ifdef HAVE_LAPACKE

#include "linearsys.h"
#include "genutil.h"
#include "explot.h"
#include <stdio.h>
#include <math.h>
#include <lapacke.h>
#include <float.h>

#define N_XS 100
#define NORMALIZATION_CONSTANT 2006.787453080206

typedef enum {
	LUP,
	MGS,
	QLESS
}Mode;

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
	
	eplot_multi_options(series, 2, &gb_settings, false);
	
	
	arr_free(&x_data_plot);
	arr_free(&y_data_plot);
}

double cond_svd(const Matrix*A) 
{
	Matrix Acp = {0};
	mat_cp(A, &Acp);
	const int ncols = A->ncols;
	const int nrows = A->nrows;
    double s[ncols], superb[ncols];
    LAPACKE_dgesvd(LAPACK_ROW_MAJOR, 'N', 'N', nrows, ncols, Acp.data, 
		ncols, s, NULL, 1, NULL, 1, superb);
	mat_free(&Acp);
    return s[0] / s[ncols-1];
}

double two_norm(const Array *a)
{
	double sum=0;
	for(int i=0; i<a->size; i++)
		sum+=pow(ARRP(a,i),2);
	return sqrt(sum);
}

double euclidean_distance(const Array *a, const Array *b)
{
	double sum=0;
	for(int i=0; i<a->size; i++)
		sum+=pow(ARRP(a,i)-ARRP(b,i),2);
	return sqrt(sum);
}

void fit_exp_sin(int n_pol, Mode mode, bool do_plot, Array *coeffs)
{
	const int n_data = 100;
	Array x_data = arr_linspace( 0, 1, n_data );
	Array y_data = arr_map( &x_data, h_norm );
	
	int *exponents = malloc(n_pol*sizeof(int));
	
	LinearModel model = linear_model_alloc(n_pol);
	for(int i=0; i<n_pol; i++)
	{
		exponents[i]=i;
		model.funcs[i].func=monomyal;
		model.funcs[i].params=exponents+i;
	}
	
	if(mode==LUP)
	{
		linst_lsqr_fit_linear_lup(&x_data, &y_data, &model, coeffs);
		Matrix A={0}, G={0};
		make_linear_model_matrix(&A, &model, &x_data);
		mat_gram(&A, &G);
		double kappa = cond_svd(&G);
		printf("Condition number of ATA: %g\n", kappa);
		/*
		Array exact={0}, atb={0};
		arr_init(&atb, y_data.size);
		mat_atb(&A, &y_data, &atb);
		linst_lsqr_fit_linear_qless(&x_data, &y_data, &model, &exact);
		printf("Two-norm difference: %g\n", euclidean_distance(coeffs, &exact));
		printf("Condition number evaluated: %g\n\n", 
			euclidean_distance(coeffs, &exact)/
				two_norm(&exact)/DBL_EPSILON);
		*/
		mat_free(&A); mat_free(&G); //arr_free(&exact); arr_free(&atb);
	}
	else if(mode==MGS)
	{
		linst_lsqr_fit_linear_qr(&x_data, &y_data, &model, coeffs);
		Matrix A={0};
		make_linear_model_matrix(&A, &model, &x_data);
		double kappa = cond_svd(&A);
		printf("Condition number of A: %g\n", kappa);
		Array exact={0};
		linst_lsqr_fit_linear_qless(&x_data, &y_data, &model, &exact);
		const double coeffs_exact_distance = euclidean_distance(coeffs, &exact);
		const double coeffs_norm = two_norm(&exact);
		printf("Two-norm difference: %g\n", coeffs_exact_distance);
		printf("Relative error  %g\n", coeffs_exact_distance/coeffs_norm);
		printf("Condition number evaluated: %g\n\n", 
			euclidean_distance(coeffs, &exact)/coeffs_norm/DBL_EPSILON);
		mat_free(&A); arr_free(&exact);
	}
	else if(mode==QLESS)
		linst_lsqr_fit_linear_qless(&x_data, &y_data, &model, coeffs);
	else
		raiseErr("Select an appropriate mode");

	

	if(do_plot)
		plot_fit( coeffs, &x_data, &y_data, &model, n_data);

	Array *arr_to_free[]={&x_data, &y_data};
	arr_free_many(arr_to_free, 2);
	free(exponents);
	linear_model_free(&model);
}

int main()
{
	Array coeffs[3]={0};
	
	
	fit_exp_sin(15, LUP,  false, coeffs);
	fit_exp_sin(15, MGS,  false, coeffs+1);
	fit_exp_sin(15, QLESS, true, coeffs+2);
	
	puts("-----------------------------------------------\n");
	
	puts(ANSI_YELLOW"Coefficients found in least squares fitting\n"ANSI_RESET);
	puts(" n |    LUP    |    MGS    |  Q-LESS");
	for(int i=0; i<15; i++)
		printf("%2d | % 9.3g | % 9.3g | % 9.3g\n", i+1, ARR(coeffs[0],i), 
			ARR(coeffs[1],i),ARR(coeffs[2],i));
	
	printf("\nQ-less n = 15 difference from 1: %g", ARR(coeffs[2],14)-1);
	
	for(int i=0; i<3; i++)
		arr_free(coeffs+i);
	
	exit(EXIT_SUCCESS);
}

#endif

#ifndef HAVE_LAPACKE
int main(){return 1;}
#endif
