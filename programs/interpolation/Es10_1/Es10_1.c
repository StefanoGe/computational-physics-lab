#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "interpolation.h"


double interpolate_cosine(int n, int k, double *cond_number)
{
	const int n_nodes = n+1;
	Array x_values = arr_linspace(-1,1,n_nodes);
	Array y_values = arr_new(n_nodes);
	Array coeffs = arr_new(n_nodes);
	
	for(int i=0; i<n_nodes; i++)
		ARR(y_values,i) = cos(k*ARR(x_values,i));

	if( !interp_polyn_vand_lup(&x_values, &y_values, &coeffs, cond_number) )
	{
		printf("for n = %d k = %d matrix was singular\n",n,k);
		return INFINITY;
	}
	
	double value_at_zero = ARR(coeffs,0);
	
	printf("For n = %d k = %d:\n",n,k);
	arr_print_inline(&coeffs, "%5.3g",true,true);
	if(cond_number) printf("cond number: %lf\n\n", *cond_number);
	
	Array *to_free[]={&x_values,&y_values,&coeffs};
	arr_free_many(to_free,3);
	return value_at_zero;
}

void test_interpolation(){
	
	const int maxk = 100;
	const int stepk = 10;
	const int nstepk = maxk/stepk;
	const int maxn = 40;
	const int stepn = 4;
	const int nstepn = maxn/stepn;

	Matrix values_at_zero = mat_new( nstepk, nstepn);
	Array cond_number = arr_new( nstepn );
	
	
//	eprint("%d %d", values_at_zero.nrows,values_at_zero.ncols);

	for( int k = stepk; k <= maxk; k+=stepk )
		for( int n = stepn; n <= maxn; n+=stepn )
		{
//			eprint("%d %d",k,n);
			MAT(values_at_zero, k/stepk-1,n/stepn-1) = interpolate_cosine( n, k, &ARR(cond_number, n/stepn-1));
		}
	mat_print_stdout(&values_at_zero,"% .3lf",true);
	puts("Logarithm of condition number of Vandermonde matrix scales with 10*n as:");
	for(int j=0; j<nstepn; j++)
	{
		printf("%5g ", log(ARR(cond_number,j)));
	}
	puts("\n");
	mat_free(&values_at_zero);
	arr_free(&cond_number);
	
}

int main()
{
	test_interpolation();

	interpolate_cosine(10,1,NULL);
	
	// Remark: condition number scales exponentially!
	
	exit(EXIT_SUCCESS);
}
