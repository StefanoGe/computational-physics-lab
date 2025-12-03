// integration.c

#include "comp_physics.h"
#include <math.h>

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

static inline double legendre_naive( double x, int degree )
{
	if(degree == 0)
		return 1;
	if(degree == 1)
		return x;
		
	const double first_addend = ( 2 * degree - 1 ) * x * legendre_naive( x, degree - 1 );
	const double second_addend = ( 1 - degree ) * legendre_naive( x, degree - 2 );
	
	return (first_addend + second_addend) / degree;
}

double legendre(double x, int degree, VectorD * pol_value_registry){return legendre_naive(x, degree);}

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

static inline double initial_point( int degree, int root_index )
{
	const double cos_arg = PI*( 4*root_index - 1 )/( 4*degree + 2);
	return cos(cos_arg);
}

double legendre_root( int degree, int root_index )
{
	Par_Func leg_par = {legendre_par, &degree, 1};
	Par_Func leg_der_par = {legendre_der_par, &degree, 1};
	
	return root_newt( &leg_par, &leg_der_par, initial_point( degree, root_index), 
				DEF_TOL, DEF_TOL, nullptr );
}


