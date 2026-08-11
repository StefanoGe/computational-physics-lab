#include "interpolation.h"
#include "explot.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

#define N_POINTS_INF_NORM 4000

double powfabs(double x, void *m){ return pown( fabs( x ), *((int*)m) ); }
/*
double one_test( Func_Ptr func, int n_nodes )
{
	BarFit barf = fit_cheb2_def( func, n_nodes );
	
	double sup_norm = inf_norm( &barf );
	
	free_vecD( &barf.points, &barf.f_values, &barf.weights );
	
	return sup_norm;
}
 
VectorD test_cheb2( Func_Ptr func , const ArrayInt * n_nodes)
{
	VectorD sup_vec = init_vecD();
	for( int i = 0; i < n_nodes->length; i ++ )
		appendD( &sup_vec, one_test(func, n_nodes-> val[i]) );
	
	return sup_vec;
}
*/

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
							int nstep)
{
	ParamFunc *fps;
	SAFE_ALLOC(fps,n_m_values);
	for(int i=0; i<n_m_values; i++)
	{
		fps[i].func = powfabs;
		fps[i].params = m_values+i;
	}
	
	const int niter = nmax_val/nstep;
	Array *inf_norms;
	SAFE_ALLOC(inf_norms,n_m_values);
	Array *log_inf_norms;
	SAFE_ALLOC(log_inf_norms,n_m_values);
	
	Array n_values = arr_linspace(nstep,nmax_val,niter);
	Array log_n_values = arr_map(&n_values,log);
	
	
	for(int i=0; i<n_m_values; i++)
	{
		inf_norms[i]=fit_one_function(fps+i,nstep,nmax_val);
		log_inf_norms[i] = arr_map(inf_norms+i,log);
	}
	
	Global gb_settings = eplot_init_global("Error study","n","||f-p||",NULL,
		"chebyshev_interp_err_study",NULL,"bottom left",NULL);
	
	SeriesSpec *series;
	SAFE_ALLOC(series,n_m_values);
	
	char label[6][20];
	
	for(int i=0; i<n_m_values; i++)
		sprintf(label[i],"m = %d", m_values[i]);
	
	for(int i = 0; i<n_m_values; i++)
		series[i] = eplot_init_series(log_n_values.data,log_inf_norms[i].data,
					niter,"lp",label[i],NULL);
	
	eplot_multi(series,n_m_values,&gb_settings);
	
	arr_free(&n_values);
	arr_free(&log_n_values);

	for(int i=0; i<n_m_values; i++)
	{
		arr_free(inf_norms+i);
		arr_free(log_inf_norms+i);
	}
	SAFE_FREE(inf_norms);
	SAFE_FREE(fps);
	SAFE_FREE(series);
}

int main()
{
	int m_values[]={1,3,5,7,9,11};
	const int n_m_values = 6;
	
	error_study_on_fabs(m_values,n_m_values,100,10);

	// Every fuction has continous every derivative up to the (m-1)th one!
	// Asymptotic behaviour should go as n^m as n goes to infinity!

	exit(EXIT_SUCCESS);
}

/*
 * OLD
 * 
	ArrayInt n_nodes = allocArrI( 10 );
	for( int i = 0; i < 10; i++ )
		n_nodes.val[i] = (i+1) * 10;
		
	VectorD n_nodesD = init_vecD();
	for( int i = 0; i < n_nodes.length; i++ )
		appendD( &n_nodesD, (double)n_nodes.val[i] );
	
	VectorD sup_norms [6];
	char * labels[6];
	
	//ArrayInt exp = allocArrI( 6 );
	for( int i = 0; i < 6; i++)
	{
		m = 1 + 2 * i;
		labels[i] = malloc( 100 * sizeof(char) );
		sprintf( labels[i], "m = %d", m );
		sup_norms[i] = test_cheb2( f1, &n_nodes );
	}
	
	FILE * gp = gp_open();
	
	char * title = malloc( sizeof(char) * 100 );
	sprintf(title, "Plot finale");
	
	fprintf(gp, "set logscale y\n");
	fprintf(gp, "set logscale x\n");
	fprintf(gp, "set key box\n");
	
	gp_term_def( gp, title );
	gp_axes_labels( gp, "n", "inf-norm error" );
	gp_set_plot(gp, 6, labels, NULL);
	for( int i = 0; i < 6; i++ )
	{
		std_print_vecD( sup_norms+i );
		gp_prt_carr( gp, n_nodesD.val, sup_norms[i].val, n_nodesD.length );
	}
	
	gp_end(gp);
	free(title);
*/

