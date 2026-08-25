#include "integration.h"
#include <math.h>

#define MAX_LEG 10'000

double intg_trap( const ParamFuncPtr f, void *params,
			double x1, double x2, int n_subint )
{
	if(n_subint<1)
		raiseErr("n_subint must be at least 1");
	
	const double interval_length = (x2-x1)/n_subint;
	
	double integral = (f(x1,params)+f(x2, params))/2;
	
	for( int i=1; i<n_subint; i++ )
		integral += f(x1 + i*interval_length,params);
	
	integral *= interval_length;
	
	return integral;
}

double intg_simp( const ParamFuncPtr f, void *params,
					double x1, double x2, int n_subint )
{
	const double interval_length = (x2 - x1) / n_subint / 2;
	
	double integral = f(x1,params) + f(x2,params);
	
	double sum_first_points = 0; 
	double sum_second_points = 0;

	for( int i = 1; i <= n_subint; i++ )
		sum_first_points += f(x1 +(2*i - 1)*interval_length, params); 
	
	for( int i = 1; i < n_subint; i++ )
		sum_second_points += f(x1 + 2*i*interval_length,params);
	
	integral += 4 * sum_first_points + 2 * sum_second_points;
	integral *= interval_length / 3;
	
	return integral;	
}


static inline void update_values( double x, int degree, double *values, 
									int stored_degrees )
{
	for( int next_n = stored_degrees; next_n<=degree; next_n++ )
	{
		const double first_addend = ( 2 * next_n - 1 ) * x * values[ next_n - 1 ];
		const double second_addend = ( 1 - next_n ) * values[ next_n - 2 ];
		values[next_n] = (first_addend + second_addend) / next_n;
	}
}

double intg_legendre( double x, int degree )
{
	static double values[MAX_LEG]={1};
	static int stored_degrees = 1;
	static double curr_x = NAN;
	
	if(degree >= MAX_LEG)
		raiseErr("Requested Legendre polynomial degree must be in [0, %d]\n"
				 "degree = %d",
		MAX_LEG-1,degree);
	
	if( x != curr_x)
	{
		values[1] = x;
		stored_degrees = 2;
		curr_x = x;
	}
		
	if( degree >= stored_degrees )
	{
		update_values( x, degree, values, stored_degrees );
		stored_degrees = degree + 1;
	}

	
	return values[degree];
}

double intg_legendre_der( double x, int degree )
{
	const double num = x*intg_legendre ( x, degree ) - intg_legendre( x, degree - 1 );
	const double den = x*x - 1;
	return degree * num/den;
}






