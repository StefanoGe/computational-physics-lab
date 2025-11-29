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
