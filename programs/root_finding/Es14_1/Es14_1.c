#include <stdio.h>
#include "rootf.h"
#include <math.h>
#include "explot.h"

#define DOMAIN_PLOT_SIZE 1000

double f1( double x, void * unused )
{
	UNUSED(unused);
	return x + 1;
}

void test1()
{	
	printf("%lf\n", root_bis( f1, NULL, -2.5, 0, 1e-10 , NULL) );
}

double f2(double x, void *unused)
{
	UNUSED(unused);
	return x*x*x - 7*x*x + 14*x -6;
}


void test2(double x1, double x2)
{
	Array debug = {0};
	arr_reserve(&debug,1);
	
	const double root = root_bis( f2, NULL, x1, x2, 0, &debug );
	printf( "Found root %.20lf\n", root );
	
	arr_print(stdout,&debug,"%.20g\n");
	
	Array log_abs_err = arr_new(debug.size);
	for(int i=0; i<debug.size; i++)
		ARR(log_abs_err,i) = -log10 ( fabs( ARR(debug,i) - root ) );
		
	Array indices = arr_linspace(1,debug.size,debug.size);

	Global gb_setting = eplot_init_global("Conergence of bisection method",
		"n","-log...",NULL,"convergence_bis",NULL,NULL,NULL);

	SeriesSpec conv = eplot_init_series(indices.data,log_abs_err.data,
		debug.size,"lp",NULL,NULL);
	
	eplot_multi(&conv,1,&gb_setting);

	arr_free(&debug);
	arr_free(&log_abs_err);
	arr_free(&indices);
}

int main()
{
	test1();
	test2(0,1);
/*
	test2(0.4, 1);
	test2(0.25, 0.66);
*/

	exit(EXIT_SUCCESS);
}

// ci dovrebbe uscire 0.5 perche' ogni volta dividiamo l'intervallo in 2
