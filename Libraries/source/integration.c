// integration.c

#include "comp_physics.h"

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

double legendre_naive( double x, int degree )
{
	if(degree == 0)
		return 1;
	if(degree == 1)
		return x;
		
	const double first_addend = ( 2 * degree - 1 ) * x * legendre_naive( x, degree - 1 );
	const double second_addend = ( 1 - degree ) * legendre_naive( x, degree - 2 );
	
	return (first_addend + second_addend) / degree;
}

double legendre(double x, int degree){return legendre_naive(x, degree);}


