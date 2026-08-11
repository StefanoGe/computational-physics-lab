#include "interpolation.h"
#include "explot.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

// capire come valutare bene questa quantità

#define N_DOMAIN_PLOT 1000

// Ricordati di studiare le proprietà della funzione
// e dimostrare che sia biunivoca (ha derivata maggiore di zero --> monotona)
// e ha come limite agli estremi di definizione +- inf

double f1(double z, void *null)
{
	(void)null;
	return 1 / ( z*z -2*z + 2 );
}
/*
void test( int num )
{
	VectorD domain = vec_range( -6, 6, DOMAIN_PLOT );
	BaricFitter baric_fitter;
	baric_fitter.points = cheb1_nodes_def( num );
	baric_fitter.weights = cheb1_weights(num);
	
	BarFit barf = bar_fir_real( &baric_fitter, f1 );
	
	VectorD y_values = init_vec_length( DOMAIN_PLOT );
	VectorD true_values = init_vec_length(DOMAIN_PLOT);
	for( int i = 0; i < DOMAIN_PLOT; i++ )
	{
		y_values.val[i] = barf_get_value_real( &barf, domain.val[i] );
		true_values.val[i] = f1(domain.val[i]);
	}
	
	VectorD y_vecs [] = { y_values, true_values };
	
	plot_mult_vecs( &domain, y_vecs, 2, NULL_INFO );
	
	free_vecD( &baric_fitter.points, &baric_fitter.weights, &barf.f_values,
				&barf.points, &barf.weights, &y_values, &true_values);
}
*/
/*
void test_inv_func ( int num)
{
	VectorD cheb1 = cheb1_nodes_def( num );
	printf( "Cheb1 points:\n" );
	std_print_vecD(&cheb1);
	double curr;
	double unit;
	
	for( int i = 0; i < num; i++)
	{
		curr = unit_to_real_line( cheb1.val[i] );
		unit = real_line_to_unit(curr);
		printf( "%10lf   %10lf  %3lf\n", curr, unit, unit/cheb1.val[i] );
	}
	
	free_vecD( &cheb1 );
}
*/

void test_realline_interp(ParamFuncPtr f,int n)
{
	ParamFunc fp = {f,0};
	BarFit barf = interp_barf_new_cheb1_points(-1,1,n);
	interp_barf_add_fvalues_realline(&barf,&fp,NULL,NULL);
	
	Array domain = arr_linspace(-6,6,N_DOMAIN_PLOT);
	Array true_values = arr_map_par(&domain,f,NULL);
	Array approx_pol = arr_map_par(&domain, interp_barf_get_value_wrap,&barf);
	
	Global gb_settings = eplot_init_global("Polynomial interpolation "
				"over the real line","x","y",NULL,
				"cheb_realline",
				NULL,NULL,"set xrange [-7:7]");
	
	SeriesSpec data = eplot_init_series(domain.data,true_values.data,
			N_DOMAIN_PLOT,"l","1/(z^2-2z+2)",NULL);
			
	SeriesSpec plot_fit = eplot_init_series(domain.data,approx_pol.data,
			N_DOMAIN_PLOT,"l dt 3","Interpolant polynomial",NULL);
	
	Array nodes = arr_new(n);
	
	for(int i=0; i<n; i++)
		ARR(nodes,i) = interp_unit_to_real_line(barf.points[i]);
			
	SeriesSpec nodes_plot = eplot_init_series(nodes.data,barf.f_values,
			n,"p","Nodes",NULL);
			
	SeriesSpec series[]={data,plot_fit,nodes_plot};
			
	eplot_multi(series,3,&gb_settings);
	
	Array *arr_to_free[]={&domain,&true_values,&approx_pol,&nodes};
	arr_free_many(arr_to_free,4);
}

int main()
{
	
	test_realline_interp(f1,30);
	
	exit(EXIT_SUCCESS);
}

