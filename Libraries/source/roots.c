#include "comp_physics.h"
#include <float.h>
#include <math.h>

#define max(a,b) (( a < b ) ? b : a)

// da sistemare tenendo conto della machine precision
static inline int sign( double x )
{
	if( x > 0 )
		return 1;
	if( x < 0 )
		return -1;
	return 0;
}

double root_bis( Par_Func fnc, double x1, double x2, double tol, VectorD * debug)
{
	const double add_tol = DBL_EPSILON * max( fabs(x1), fabs(x2) );
	int sign1 = sign( evaluate( &fnc, x1 ) );
	int sign_curr;
	double curr;
	
	while( fabs( x1 - x2 ) > tol + add_tol )
	{
		curr = (x1 + x2) / 2;
		sign_curr = sign( evaluate( &fnc, curr) );
		if( sign1 == sign_curr )
		{
			x1 = curr;
			sign1 = sign_curr;
		} else if( sign_curr == 0 )
			return curr;
		else
			x2 = curr;
		if(debug)
			appendD( debug, curr );
	}
	return curr;
}

