#include "comp_physics.h"
#include <float.h>
#include <math.h>
#include <stdbool.h>

#define max(a,b) (( a < b ) ? b : a)

#define DEFAULT_TOLERANCE 100*DBL_EPSILON
#define MAX_ITERATION 10000

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



double root_newt( const Par_Func * fnc, const Par_Func * derivative, 
					double start, double xtol, double ftol, VectorD * debug )
{
	if( xtol == DEF_TOL )
		xtol = DEFAULT_TOLERANCE;
		
	if( ftol == DEF_TOL )
		ftol = DEFAULT_TOLERANCE;
	
	bool found_root = false;
	
	double prev_x = 0;
	double curr_x = start;
	
	if(debug)
		if(debug->val)
		{
			free_vecD( debug );
			*debug = init_vecD();
			appendD( debug, start );
		}
	
	for( int i = 0; i < MAX_ITERATION && !found_root; i++ )
	{
		prev_x = curr_x;
		curr_x -= evaluate( fnc, curr_x ) / evaluate( derivative, curr_x );
		if(debug)
			appendD( debug, curr_x );
		found_root = fabs( prev_x - curr_x ) <= xtol 
					&& fabs( evaluate( fnc, curr_x ) ) <= ftol;
	}
	
	if( !found_root )
		raiseErr( "Loop did not stop after %d iterations", MAX_ITERATION );
	
	return curr_x;
}





