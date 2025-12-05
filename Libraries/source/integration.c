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
// scales bad with degree
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

static inline void leg_smart_reset( int * length, double * curr_x, double new_x  )
{
	*length = 0;
	*curr_x = new_x;
}

static inline double legendre_smart( double x, int degree );

static inline void update_values( double x, int degree, double * values )
{
	if(degree == 0)
		values[degree] = 1;
	
	else if(degree == 1)
		values[degree] = x;
	
	else
	{
		const double first_addend = ( 2 * degree - 1 ) * x * legendre_smart( x, degree - 1 );
		const double second_addend = ( 1 - degree ) * legendre_smart( x, degree - 2 );
		
		values[degree] = (first_addend + second_addend) / degree;
	}
}

// scales as O(degree)
static inline double legendre_smart( double x, int degree )
{
	static double values[MAX_LEG];
	static int length = 0;
	static double curr_x = NAN;
	
	if( x != curr_x)
		leg_smart_reset( &length, &curr_x, x );
		
	if( degree >= length )
		update_values( x, degree, values );
	
	return values[degree];
}


double legendre(double x, int degree){return legendre_naive(x, degree);}

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
				DBL_EPSILON*10, DBL_EPSILON*10, nullptr );
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
		eprint( "root ok" );
		const double curr_weight = gauss_legendre_weight( curr_root, order );
		eprint("weight ok");
		const double fnc_value = evaluate( fnc, curr_root );
		eprint("fnc_value ok");
		sum += curr_weight * fnc_value;
	}
	return sum;
}






