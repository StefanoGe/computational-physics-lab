#include "comp_physics.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define N_DOMAIN 1000

double cos_sin(double x) { return cosh(sin(x)); }

void test_cheb( Func_Ptr func, int n_nodes, double x1, double x2 )
{
	VectorD domain = vec_range( x1, x2, N_DOMAIN );
	BarFit barf = fit_cheb2( func, n_nodes, x1, x2 );
	VectorD y_values = init_vecD();
	VectorD true_values = init_vecD();
	for( int i = 0; i < N_DOMAIN; i++ )
	{
		appendD( &y_values, barf_get_value( &barf, domain.val[i] ) );
		appendD( &true_values, func(domain.val[i]) );
	}
	
	VectorD ys[] = { y_values, true_values };
	
	PlotInfo plot_info = plot_info_init(2);
	plot_info.title = "'Interpolation of cosh(sin(x))'";
	plot_info.labels[0] = "Fitted function";
	plot_info.labels[1] = "cosh(sin(x))";
	
	plot_mult_vecs( &domain, ys, 2, plot_info );
	
	free_vecD( &y_values, &true_values, &domain, &barf.f_values, &barf.points,
	&barf.weights );
}

int main()
{
	
	test_cheb( cos_sin, 40, 0, 2 * PI );
	
	exit(EXIT_SUCCESS);
}

