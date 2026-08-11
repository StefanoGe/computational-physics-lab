#include "rootf.h"
#include "float.h"
#include <math.h>

static inline int sign( double x )
{
	if( x > 0 )
		return 1;
	if( x < 0 )
		return -1;
	return 0;
}

double root_bis( ParamFuncPtr fnc, void *params, 
			double x1, double x2, double tol, Array *debug)
{
	const double add_tol = DBL_EPSILON * MAX( fabs(x1), fabs(x2) );
	int sign1 = sign( fnc(x1,params) );
	int sign_curr = 0;
	double curr = 0;
	
	while( fabs( x1 - x2 ) > tol + add_tol )
	{
		curr = (x1 + x2) / 2;
		sign_curr = sign( fnc(curr,params) );
		if( sign1 == sign_curr )
		{
			x1 = curr;
			sign1 = sign_curr;
		} else if( sign_curr == 0 )
			return curr;
		else
			x2 = curr;
		if(debug)
			arr_append( debug, curr );
	}
	return curr;
}
