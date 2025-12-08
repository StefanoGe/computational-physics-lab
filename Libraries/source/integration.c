// integration.c
#include "comp_physics.h"
#include <math.h>
#include <float.h>

#define MAX_LEG 100

double int_trap( const Par_Func * fnc, double x1, double x2, int n_subint )
{
	const double interval_length = (x2 - x1) / n_subint;
	
	double integral = ( evaluate(fnc, x1) + evaluate(fnc, x2) )/2;
	
	for( int i = 1; i < n_subint; i++ )
		integral += evaluate( fnc, x1 + i * interval_length );
	
	integral *= interval_length;
	
	return integral;
}

double int_simp( const Par_Func * fnc, double x1, double x2, int n_subint )
{
	const double interval_length = (x2 - x1) / n_subint / 2;
	
	double integral = evaluate(fnc, x1) + evaluate(fnc, x2);
	
	double sum_first_points = 0; 
	double sum_second_points = 0;

	for( int i = 1; i <= n_subint; i++ )
		sum_first_points += evaluate( fnc, x1 + ( 2 * i - 1 ) * interval_length ); 
	
	for( int i = 1; i < n_subint; i++ )
		sum_second_points += evaluate( fnc, x1 + 2 * i * interval_length );

	
	integral += 4 * sum_first_points + 2 * sum_second_points;
	integral *= interval_length / 3;
	
	return integral;	
}

static inline void update_values( double x, int degree, double *values, 
									int length )
{
	for( int next_n = length; next_n<=degree; next_n++ )
	{
		const double first_addend = ( 2 * next_n - 1 ) * x * values[ next_n - 1 ];
		const double second_addend = ( 1 - next_n ) * values[ next_n - 2 ];
		values[next_n] = (first_addend + second_addend) / next_n;
	}
}

double legendre( double x, int degree )
{
	static double values[MAX_LEG]={1};
	static int length = 1;
	static double curr_x = NAN;
	
	if( x != curr_x)
	{
		values[1] = x;
		length = 2;
		curr_x = x;
	}
		
	if( degree >= length )
	{
		update_values( x, degree, values, length );
		length = degree + 1;
	}

	
	return values[degree];
}


double legendre_par(double x, void *degree) {
    int d = *(int *)degree;
    return legendre(x, d);
}

double legendre_der( double x, int degree )
{
	const double num = x*legendre ( x, degree ) - legendre( x, degree - 1 );
	const double den = x*x - 1;
	return degree * num/den;
}

double legendre_der_par(double x, void *degree) {
    int d = *(int *)degree;
    return legendre_der(x, d);
}

double cheb_roots( int degree, int root_index )
{
	const double cos_arg = PI*( 4*root_index - 1 )/( 4*degree + 2);
	return cos(cos_arg);
}

double legendre_root( int degree, int root_index )
{
	eprint("asked for degree %d and index %d", degree, root_index);
	Par_Func leg_par = {legendre_par, &degree, 1};
	Par_Func leg_der_par = {legendre_der_par, &degree, 1};
	
	return root_newt( &leg_par, &leg_der_par, cheb_roots( degree, root_index), 
				10*DBL_EPSILON, 10*DBL_EPSILON, nullptr );
}

static inline double gauss_legendre_weight( double x, int n )
{
	const double num = 1 - x * x;
	const int den1 = ( n + 1 ) * ( n + 1 );
	const double den2 = pown( legendre( x, n + 1 ), 2 );
	return 2 * num / den1 / den2;
}

double int_gauss_legendre( Par_Func * fnc, int order )
{
	double sum = 0;
	double curr_root = NAN;
	for( int i = 1; i <= order; i++)
	{
		
		curr_root = legendre_root( order, i );
		const double curr_weight = gauss_legendre_weight( curr_root, order );
		const double fnc_value = evaluate( fnc, curr_root );
		sum += curr_weight * fnc_value;
	}
	return sum;
}






