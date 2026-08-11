#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "interpolation.h"
#include "explot.h"

#define PLOT_DOMAIN_LENGTH 100
#define N_POINTS_INF_NORM 4000

double f1(double x, void*null)
{
	(void)null;
	return 1.0/(25*x*x+1);
}

double f2(double x, void*null)
{
	(void)null;
	return tanh(5*x+2);
}

double f3(double x, void*null)
{
	(void)null;
	return cosh(sin(x));
}

double f4(double x, void*null)
{
	(void)null;
	return sin(cosh(x));
}

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
							- barf->func( curr_point,NULL ) );
		if( max_diff < curr_diff )
			max_diff = curr_diff;
	}
	return max_diff;
}

Array fit_one_function(ParamFuncPtr f, int step, int max_val)
{
	const int niter = max_val / step;
	Array inf_norms = arr_new(niter);
	ParamFunc fp={f,NULL};
	for(int n=step,i=0; n<=max_val; n+=step, i++)
	{
		BarFit barf = interp_barf_new_cheb2_points(-1,1,n);
		interp_barf_add_fvalues(&barf,&fp);
		ARR(inf_norms,i) = inf_norm(&barf);
		interp_barf_free(&barf);
	}
	
//	arr_print_inline(&inf_norms, "%3g ",true,true);
	
	return inf_norms;
}

void collect_and_plot(ParamFuncPtr *fs, int count,int step, int max_val,
		char **f_names)
{
	const int niter = max_val/step;
	Array *inf_norm_errs;
	SAFE_ALLOC(inf_norm_errs,count);
	Array n_values = arr_linspace(4,max_val,niter);
	
	for(int i = 0; i<count; i++)
		inf_norm_errs[i] = fit_one_function(fs[i],step,max_val);
	
	Global gb_settings = eplot_init_global("Error study","n","||f-p||",NULL,
		"chebyshev_interp_err_study",NULL,"bottom left",NULL);
	
	SeriesSpec *series;
	SAFE_ALLOC(series,count);
	
	Array *log_inf_norms;
	SAFE_ALLOC(log_inf_norms,count);
	
	for(int i = 0; i<count; i++)
	{
		log_inf_norms[i]=arr_map(inf_norm_errs+i,log);
		series[i] = eplot_init_series(n_values.data,log_inf_norms[i].data,
					niter,"lp",f_names[i],NULL);
	}
	
	eplot_multi(series,count,&gb_settings);
	
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
	ParamFuncPtr fs[] = {f1,f2,f3,f4};
	const int count = 4;
	collect_and_plot(fs,count,4,60,f_names);
	
	// ad an estimate after you have refractored least squares interpolation
	// and add same for equidistant nodes

	exit(EXIT_SUCCESS);
}
