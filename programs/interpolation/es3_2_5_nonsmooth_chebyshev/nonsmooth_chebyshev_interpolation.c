#include "interpolation.h"
#include "explot.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "linearsys.h"

#define N_POINTS_INF_NORM 4000

double powfabs(double x, void *m){ return pown( fabs( x ), *((int*)m) ); }

double one(double x, void*unused){UNUSED(x); UNUSED(unused); return 1;}

double id(double x, void*unused){UNUSED(unused); return x;}

double line(double x, void *params){
	return ((double*)params)[0]+((double*)params)[1]*x;}

double inf_norm( const BarFit * barf)
{
	double max_diff = 0;
	const double length = 2;
	double curr_diff;
	double curr_point;
	for( int i = 0; i < N_POINTS_INF_NORM; i ++ )
	{
		curr_point = -1 + (length * i) /N_POINTS_INF_NORM;
		curr_diff = fabs( interp_barf_get_value( barf, curr_point ) 
							- barf->func( curr_point,barf->params ) );
		if( max_diff < curr_diff )
			max_diff = curr_diff;
	}
	return max_diff;
}

Array fit_one_function(const ParamFunc *f, int step, int max_val)
{
	const int niter = max_val / step;
	Array inf_norms = arr_new(niter);
	for(int n=step,i=0; n<=max_val; n+=step, i++)
	{
		BarFit barf = interp_barf_new_cheb2_points(-1,1,n);
		interp_barf_add_fvalues(&barf,f);
		ARR(inf_norms,i) = inf_norm(&barf);
		interp_barf_free(&barf);
	}
	
//	arr_print_inline(&inf_norms, "%3g ",true,true);
	
	return inf_norms;
}

void error_study_on_fabs(int m_values[], int n_m_values, int nmax_val,
							int nstep, const int points_in_fit[])
{
	ParamFunc *fps;
	SAFE_ALLOC(fps,n_m_values);
	for(int i=0; i<n_m_values; i++)
	{
		fps[i].func = powfabs;
		fps[i].params = m_values+i;
	}
	
	const int niter = nmax_val/nstep;
	Array *inf_norm_errs;
	SAFE_ALLOC(inf_norm_errs,n_m_values);
	Array *log_inf_norms;
	SAFE_ALLOC(log_inf_norms,n_m_values);
	
	Array n_values = arr_linspace(nstep,nmax_val,niter);
	Array log_n_values = arr_map(&n_values, log10);
	
	Array *fit_coeffs;
	SAFE_ALLOC(fit_coeffs, n_m_values);
	
	Array *fitted_err_values;
	SAFE_ALLOC(fitted_err_values, n_m_values);
	
	ParamFuncPtr funcs_to_fit[] = {one, id};
	
	
	SeriesSpec *series;
	SAFE_ALLOC(series, n_m_values*2);
	
	char (*label)[20];
	SAFE_ALLOC(label, n_m_values);
	
	for(int i=0; i<n_m_values; i++)
		sprintf(label[i],"m = %d", m_values[i]);
	
	char (*linestyle)[50];
	SAFE_ALLOC(linestyle,n_m_values);
	
	for(int i = 0; i<n_m_values; i++)
	{
		snprintf(linestyle[i], 50, "lp ls %d", i+1);
		
		inf_norm_errs[i]=fit_one_function(fps+i,nstep,nmax_val);
		log_inf_norms[i]=arr_map(inf_norm_errs+i,log10);
		series[2*i] = eplot_init_series(log_n_values.data,log_inf_norms[i].data,
					niter, linestyle[i], label[i], NULL);
		arr_null(fit_coeffs+i);
		
		const int n_points_to_fit = points_in_fit[i*2+1]+1-points_in_fit[i*2];
		
		Array temp_inf_norm = {log_inf_norms[i].data + points_in_fit[2*i], 
				n_points_to_fit, n_points_to_fit, false};
				
		Array temp_n_values = {log_n_values.data + points_in_fit[2*i], 
				n_points_to_fit, n_points_to_fit, false};
		
		linst_lsqr_fit_linear_qr_wrap(&temp_n_values, &temp_inf_norm, funcs_to_fit, 
			NULL, 2, fit_coeffs+i);
		
		printf(ANSI_GREEN"For the function "ANSI_RED"%-13s "ANSI_GREEN"the "
			"constant K is: "ANSI_YELLOW"%.3g"ANSI_RESET"\n",
				"|x|^m_values[i]", -ARRP(fit_coeffs+i,1));
		
		fitted_err_values[i] = arr_map_par(&temp_n_values, line, 
			(void*)(fit_coeffs+i)->data);

		series[i*2+1] = eplot_init_series(temp_n_values.data, 
					fitted_err_values[i].data, n_points_to_fit, "l dt 2", NULL,NULL);
	}
	
	Global gb_settings = eplot_init_global("Error study","log10(n)",
		"log10(||f-p||)",NULL,
		"chebyshev_interp_err_study",NULL,"bottom left",NULL);

	
	eplot_multi_options(series, n_m_values*2, &gb_settings, false);
	
	arr_free(&n_values);
	arr_free(&log_n_values);

	for(int i=0; i<n_m_values; i++)
	{
		arr_free(inf_norm_errs+i);
		arr_free(log_inf_norms+i);
	}
	SAFE_FREE(inf_norm_errs);
	SAFE_FREE(fps);
	SAFE_FREE(series);
}

int main()
{
	int m_values[]={1,3,5,7,9,11};
	const int n_m_values = 6;
	const int points_in_fit[] = {0, 9, 0, 9, 1, 9, 1, 9, 1, 9, 1, 8};
	
	error_study_on_fabs(m_values, n_m_values, 100, 10, points_in_fit);

	// Every fuction has continous every derivative up to the (m-1)th one!
	// Asymptotic behaviour should go as n^m as n goes to infinity!

	exit(EXIT_SUCCESS);
}

