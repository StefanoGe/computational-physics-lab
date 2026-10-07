#include "integration.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "explot.h"

#define MAX_DEGREE_TEST 5
#define DOMAIN_PLOT_SIZE 100

const double leg_roots[3][4] = {
	{ 0.5773502691896258, -0.5773502691896258 },
	{ 0.7745966692414834, 0, -0.7745966692414834 },
	{ 0.8611363115940526, 0.3399810435848563,
	 -0.3399810435848563, -0.8611363115940526 }
};

const double leg_weights[3][4] = {
	{1,1},
	{0.5555555555555555555555555555,0.8888888888888888888888,
		0.55555555555555555555555555555},
	{0.347854845137453857373,0.652145154862546142627,
			0.652145154862546142627,0.347854845137453857373,}
};

/*

void plot_polynomials()
{
	VectorD xaxis = vec_range( -1, 1, 100);
	MatrixDouble yvalues = allocMatD( 2, 100 );
	for( int i = 0; i < 100; i++)
	{
		double x = xaxis.val[i];
		yvalues.val[0][i] = legendre( x, 2 );
		yvalues.val[1][i] = ( 3*x*x - 1 ) / 2;
	}
	
	tmplot_carrs( "plot_pol", xaxis.val, yvalues.val, 100, 2 );
	
	freeMatD(yvalues);
	free_vecD(&xaxis);
}


void plot_polynomials24()
{
	int n_points = 1000;
	
	VectorD xaxis = vec_range( -1, 1, n_points);
	MatrixDouble yvalues = allocMatD( 1, n_points );
	for( int i = 0; i < n_points; i++)
	{
		double x = xaxis.val[i];
		yvalues.val[0][i] = legendre( x, 24 );
	}
	
	tmplot_carrs( "plot_pol", xaxis.val, yvalues.val, n_points, 1 );
	
	freeMatD(yvalues);
	free_vecD(&xaxis);
}
*/

void test_roots(){
	
	Array debug={0};
	
	for( int deg = 2; deg <= MAX_DEGREE_TEST; deg++ )
	{
		printf( ANSI_GREEN"Roots and weights for Gauss quadrature rule "
				"of Legendre polynomial of degree "ANSI_RED"%d\n"
				ANSI_YELLOW"Root index | Root value"
				" | Weight"
				/*" | Newt steps"*/, deg);
		if( 2 <=deg && deg <=4 )
			printf(" | Diff with true value");
		putchar('\n');
		
		if( 2 <=deg && deg <=4 )
		{		
			const int needed_space=33;
			for(int i=0; i<needed_space; i++)
				putchar(' ');
			printf("| Roots    | Weights");
		}
		
		putchar('\n');
		printf(ANSI_RESET);
		for( int root = 1; root <= deg; root++ )
		{
			double curr_root = intg_legendre_root(deg, root, &debug);
			double curr_weight = intg_gauss_legendre_weight(curr_root,deg);
			printf( "%10d | % .7lf"
					" | %.4lf"
			/*"| %-10d"*/, root, curr_root, curr_weight
			/*, debug.size*/ );
		if( 2 <=deg && deg <=4 )
			printf(" | %.2e | %.2e", fabs(curr_root-leg_roots[deg-2][root-1]),
				fabs(curr_weight-leg_weights[deg-2][root-1]));
		putchar('\n');
		}
		puts("\n-------------------------------------\n");
	}
	
	arr_free(&debug);
}

void plot_polynomials(int max_deg)
{
	Array domain = arr_linspace(-1,1,DOMAIN_PLOT_SIZE);
	Array *y_leg_values;
	SAFE_ALLOC(y_leg_values,max_deg+1);
	for(int i=0; i<max_deg+1; i++)
		arr_null(y_leg_values+i);
	SeriesSpec *series;
	SAFE_ALLOC(series,max_deg+1);
	char *labels;
	SAFE_ALLOC(labels,(max_deg+1)*100);
	
	for(int deg=0; deg<=max_deg; deg++)
	{
		snprintf(labels+100*deg,100,"P_{%d}(x)",deg);
		y_leg_values[deg]=arr_map_par(&domain,intg_legendre_wrap,(void*)&deg);
		series[deg]=eplot_init_series(domain.data,(y_leg_values+deg)->data,
			domain.size,"l",labels+100*deg,NULL);
	}
	
	Global gb_settings = eplot_init_global("Legendre polynomials",
		"x","y",NULL,"leg_pol",NULL,"bottom right",
		"set xrange [-1.1:1.3]\nset yrange [-1.1:1.1]\n");
	
	eplot_multi(series,max_deg+1,&gb_settings);
	
	SAFE_FREE(series);
	arr_free(&domain);
	for(int i=0; i<max_deg+1; i++)
		arr_free(y_leg_values+i);
	SAFE_FREE(y_leg_values);
	SAFE_FREE(labels);
}

// Remark: it was necessary to set 1000 DBL epsilon in ftol. A safer 
// implementation would set them higher. 
// Remark: legendre memoization

int main()
{
//	eprint("setup ok");
//	plot_polynomials24();
//	eprint("plot ok");
	test_roots();
	plot_polynomials(5);

	
	exit(EXIT_SUCCESS);
}

