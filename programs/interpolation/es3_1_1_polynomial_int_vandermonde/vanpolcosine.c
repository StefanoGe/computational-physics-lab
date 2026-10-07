#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "interpolation.h"
#include "explot.h"

double interpolate_cosine(int n, int k, double *cond_number)
{
	const int n_nodes = n+1;
	Array x_values = arr_linspace(-1,1,n_nodes);
	Array y_values = arr_new(n_nodes);
	Array coeffs = arr_new(n_nodes);
	
	for(int i=0; i<n_nodes; i++)
		ARR(y_values,i) = cos(k*ARR(x_values,i));

	if( interp_polyn_vand_lup(&x_values, &y_values, &coeffs, cond_number) == SINGULAR)
	{
		printf("for n = %d k = %d matrix was singular\n",n,k);
		return INFINITY;
	}
	
	double value_at_zero = ARR(coeffs,0);
	
//	printf("For n = %d k = %d:\n",n,k);
//	arr_print_inline(&coeffs, "%5.3g",true,true);
//	if(cond_number) printf("cond number: %lf\n\n", *cond_number);
	
	Array *to_free[]={&x_values,&y_values,&coeffs};
	arr_free_many(to_free,3);
	return value_at_zero;
}

void plot_cond_number(const Array *cond_number)
{
	Array n_values = arr_linspace(5,41,10);
	
	SeriesSpec series = eplot_init_series(n_values.data, cond_number->data,
		cond_number->size, "lp", NULL, NULL);
	
	Global gb_settings = eplot_init_global(
		"Condition number of Vandermonde matrix", "System size", 
		"Condition number", "y", NULL, "vandermonde_cond_number", "","");
	
	eplot_multi(&series, 1, &gb_settings);
	
	arr_free(&n_values);
}

void test_interpolation(){
	
	const int maxk = 100;
	const int stepk = 10;
	const int nstepk = maxk/stepk;
	const int maxn = 40;
	const int stepn = 4;
	const int nstepn = maxn/stepn;

	Matrix values_at_zero = mat_new( nstepk, nstepn);
	Matrix one_matrix = mat_new( nstepk, nstepn);
	mat_setv(&one_matrix,1);
	Array cond_number = arr_new( nstepn );
	
	
//	eprint("%d %d", values_at_zero.nrows,values_at_zero.ncols);

	for( int k = stepk; k <= maxk; k+=stepk )
		for( int n = stepn; n <= maxn; n+=stepn )
		{
//			eprint("%d %d",k,n);
			MAT(values_at_zero, k/stepk-1,n/stepn-1) = interpolate_cosine( n, k, &ARR(cond_number, n/stepn-1));
		}
	
	puts("-----------------------------------------------------------------------------------"
		 "------------------\n");
	
	mat_axpy(&values_at_zero, -1, &one_matrix);
	puts(ANSI_GREEN"Values of fitting polynomials in x = 0\n"
		 "Function to fit: "ANSI_RED"cos(kx)\n"
		 ANSI_YELLOW"Expected value: 1\n"
		 "Below are reported the differences between the values found and the "
		 "expected value.\n"
		 "n is the number of equidistant points in [-1:1] used in the interpolation\n"ANSI_RESET);
	
	fputs("   n  ", stdout);
	
	for(int i=0; i<10; i++)
		printf("%8d ",(i+1)*4+1);
		
	puts("\nk\n");
	
	for(int row=0; row<values_at_zero.nrows; row++)
	{
		printf("%3d   ", 10*(row+1) );
		for(int col=0; col<values_at_zero.ncols; col++)
		{
			printf("% 8.1g ", MAT(values_at_zero,row,col));
		}
		putchar('\n');
	} 
	
	putchar('\n');
	
	puts("-----------------------------------------------------------------------------------"
		 "------------------\n");


	plot_cond_number(&cond_number);

	mat_free(&values_at_zero);
	arr_free(&cond_number);
	mat_free(&one_matrix);
}

int main()
{
	test_interpolation();

	interpolate_cosine(10,1,NULL);
	
	// Remark: condition number scales exponentially!
	
	exit(EXIT_SUCCESS);
}
