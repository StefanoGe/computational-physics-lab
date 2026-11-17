#include <stdio.h>
#include "genutil.h"
#include <math.h>
#include "array.h"
#include "explot.h"
#include <float.h>

double f1(double x)
{
	return (exp(x)-1)/x;
}

double f1_maclaurin(double x, int n)
{
	double sum=1;
	for(int i=0; i<n; i++)
	{
		sum *= x/(n+1-i);
		sum += 1;
	}
	return sum;
}

int main()
{
	puts("x     | Naive evaluation   | Maclaurin series   | Difference");
	const int starting_exponent =  3;
	const int last_exponent     = 16;
	const int n_calculations    = last_exponent - starting_exponent+1;
	
	Array x_values = arr_new(n_calculations);
	Array differences = arr_new(n_calculations);
	
	for(int i=starting_exponent; i<=last_exponent; i++)
	{
		const double x = pow(10,-i);
		ARR(x_values,i-starting_exponent)=x;
		const double f1_value = f1(x);
		const double f1_maclaurin_value = f1_maclaurin(x, 4);
		ARR(differences, i-starting_exponent) = fabs(f1_value-f1_maclaurin_value)/DBL_EPSILON;
		printf("%.0e | %.16lf | %.16lf | %+e\n", x,
			f1_value, f1_maclaurin_value,
			f1_value-f1_maclaurin_value);
	}
	
	Global gb_settings = eplot_init_global("Error of (e^x-1)/",
		"x", "error/ϵ_{mach}",
		"xy","another_exponential_error","another_exponential_error",
		"b r",  "set xrange[0.01:1e-17]");
	
	SeriesSpec series = eplot_init_series(x_values.data, differences.data,
		n_calculations, "lp", NULL,NULL);
	
	eplot_multi_options(&series, 1, &gb_settings, false);
	
	arr_free(&x_values);
	arr_free(&differences);
	
	return 0;
}
