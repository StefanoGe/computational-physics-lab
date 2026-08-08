#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "interpolation.h"

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

void paperino(ParamFuncPtr f, int step, int max_val)
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
	
	arr_print_inline(&inf_norms, "%3g ",true,true);
	
	Array n_values = arr_linspace(0,niter-1,niter);
	
	
	
	arr_free(&inf_norms);
	arr_free(&n_values);
}



int main()
{
	
	paperino(f1,4,60);
	

	exit(EXIT_SUCCESS);
}
