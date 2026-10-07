#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "interpolation.h"
#include "explot.h"
#include "linearsys.h"

#define PLOT_DOMAIN_LENGTH 100
#define N_POINTS_INF_NORM 4000

typedef enum {
	CHEBYSHEV,
	EQUIDISTANT
}AnalysisMode;

double f1(double x, void*null){
	(void)null;
	return 1.0/(25*x*x+1);}

double f2(double x, void*null){
	(void)null;
	return tanh(5*x+2);}

double f3(double x, void*null){
	(void)null;
	return cosh(sin(x));}

double f4(double x, void*null){
	(void)null;
	return sin(cosh(x));}

double one(double x, void*unused){UNUSED(x); UNUSED(unused); return 1;}

double id(double x, void*unused){UNUSED(unused); return x;}

double line(double x, void *params){
	return ((double*)params)[0]+((double*)params)[1]*x;}

double inf_norm( const BarFit * barf )
{
	double max_diff = 0;
	const double length = 2;
	double curr_diff;
	double curr_point;
	for( int i = 0; i < N_POINTS_INF_NORM; i ++ )
	{
		curr_point = -1 + (length * i) /N_POINTS_INF_NORM;
		curr_diff = fabs( interp_barf_get_value( barf, curr_point ) 
							- barf->func( curr_point, NULL ) );
		if( max_diff < curr_diff )
			max_diff = curr_diff;
	}
	return max_diff;
}

Array fit_one_function(ParamFuncPtr f, int step, int max_val, AnalysisMode mode)
{
	const int niter = max_val / step;
	Array inf_norms = arr_new(niter);
	ParamFunc fp={f,NULL};
	for(int n=step,i=0; n<=max_val; n+=step, i++)
	{
		BarFit barf;
		if(mode == CHEBYSHEV)
			barf = interp_barf_new_cheb2_points(-1,1,n);
		else if(mode == EQUIDISTANT)
			barf = interp_barf_new_eq_points(-1,1,n);
		interp_barf_add_fvalues(&barf,&fp);
		ARR(inf_norms,i) = inf_norm(&barf);
		interp_barf_free(&barf);
	}
	
//	arr_print_inline(&inf_norms, "%3g ",true,true);
	
	return inf_norms;
}

void collect_and_plot_chebyshev(ParamFuncPtr *fs, int count,int step, int max_val,
		char **f_names, int points_in_fit[])
{
	const int niter = max_val/step;
	Array *inf_norm_errs;
	SAFE_ALLOC(inf_norm_errs,count);
	Array n_values = arr_linspace(4,max_val,niter);
	
	for(int i = 0; i<count; i++)
		inf_norm_errs[i] = fit_one_function(fs[i], step, max_val, CHEBYSHEV);
	
	Global gb_settings = eplot_init_global("Chebyshev nodes fit error","n",
		"||f-p|| (log scale)", NULL,
		"chebyshev_interp_err_study",
		"chebyshev_interp_err_study","bottom left", "set yrange [-17:1]");
	
	SeriesSpec *series;
	SAFE_ALLOC(series,count*2);
	
	Array *log_inf_norms;
	SAFE_ALLOC(log_inf_norms,count);
	
	Array *fit_coeffs;
	SAFE_ALLOC(fit_coeffs, count);
	
	Array *fitted_err_values;
	SAFE_ALLOC(fitted_err_values, count);
	
	ParamFuncPtr funcs_to_fit[] = {one, id};
	
	for(int i = 0; i<count; i++)
	{
		log_inf_norms[i]=arr_map(inf_norm_errs+i,log10);
		series[i*2] = eplot_init_series(n_values.data,log_inf_norms[i].data,
					niter,"lp",f_names[i],NULL);
		arr_null(fit_coeffs+i);
		
		Array temp_inf_norm = {log_inf_norms[i].data, points_in_fit[i], 
				points_in_fit[i], false};
				
		Array temp_n_values = {n_values.data, points_in_fit[i], 
				points_in_fit[i], false};
		
		linst_lsqr_fit_linear_qr_wrap(&temp_n_values, &temp_inf_norm, funcs_to_fit, 
			NULL, 2, fit_coeffs+i);
		
		printf(ANSI_GREEN"For the function "ANSI_RED"%-13s "ANSI_GREEN"the "
			"constant K is: "ANSI_YELLOW"%.3g"ANSI_RESET"\n",
				f_names[i], pow(10,-ARRP(fit_coeffs+i,1)));
		
		fitted_err_values[i] = arr_map_par(&n_values, line, 
			(void*)(fit_coeffs+i)->data);

		series[i*2+1] = eplot_init_series(n_values.data, 
					fitted_err_values[i].data, niter, "l dt 2", NULL,NULL);
	}
	
	eplot_multi_options(series,count*2, &gb_settings, false);
	
	for(int i=0; i<count; i++)
	{
		arr_free(log_inf_norms+i);
		arr_free(inf_norm_errs+i);
		arr_free(fit_coeffs+i);
		arr_free(fitted_err_values+i);
	}
	SAFE_FREE(inf_norm_errs);
	arr_free(&n_values);
	SAFE_FREE(series);
	SAFE_FREE(fit_coeffs);
	SAFE_FREE(fitted_err_values);
}

void collect_and_plot_equidistant(ParamFuncPtr *fs, int count,int step, int max_val,
		char **f_names)
{
	const int niter = max_val/step;
	Array *inf_norm_errs;
	SAFE_ALLOC(inf_norm_errs,count);
	Array n_values = arr_linspace(4,max_val,niter);
	
	for(int i = 0; i<count; i++)
		inf_norm_errs[i] = fit_one_function(fs[i], step, max_val, EQUIDISTANT);
	
	Global gb_settings = eplot_init_global("Equidistant nodes fit error","n",
		"||f-p|| (log scale)", NULL,
		"equidistant_interp_err_study",
		"equidistant_interp_err_study","bottom left", "set yrange [-17:10]");
	
	SeriesSpec *series;
	SAFE_ALLOC(series,count);
	
	Array *log_inf_norms;
	SAFE_ALLOC(log_inf_norms,count);
	
	for(int i = 0; i<count; i++)
	{
		log_inf_norms[i]=arr_map(inf_norm_errs+i,log10);
		series[i] = eplot_init_series(n_values.data,log_inf_norms[i].data,
					niter,"lp",f_names[i],NULL);
	}
	
	eplot_multi_options(series,count, &gb_settings, false);
	
	for(int i=0; i<count; i++)
	{
		arr_free(log_inf_norms+i);
		arr_free(inf_norm_errs+i);
	}
	SAFE_FREE(inf_norm_errs);
	arr_free(&n_values);
	SAFE_FREE(series);
}

int main()
{
	char *f_names[]={"1/(25x^2 + 1)","tanh(5x+2)","cosh(sin x)","sin(cosh x)"};
	ParamFuncPtr fs[] = {f1, f2, f3, f4};
	const int count = 4;
	int points_in_fit[]={15,15,5,5};
	collect_and_plot_chebyshev(fs, count, 4, 60, f_names, points_in_fit);
	collect_and_plot_equidistant(fs, count, 4, 60, f_names);


	// and add same for equidistant nodes

	exit(EXIT_SUCCESS);
}
