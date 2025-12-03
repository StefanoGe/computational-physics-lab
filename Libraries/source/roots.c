#include "comp_physics.h"
#include <float.h>
#include <math.h>
#include <stdbool.h>

#define max(a,b) (( a < b ) ? b : a)

#define DEFAULT_TOLERANCE DBL_EPSILON
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

static inline void set_default_tolerance_if_requested( double * xtol, double * ftol)
{
	if( *xtol == DEF_TOL )
		*xtol = DEFAULT_TOLERANCE;
		
	if( *ftol == DEF_TOL )
		*ftol = DEFAULT_TOLERANCE;
}

static inline void init_debug( VectorD * debug )
{
	if(debug->val)
	{
		free_vecD( debug );
		*debug = init_vecD();
	}
}

static inline void newt_next_step( const Par_Func * fnc, const Par_Func * derivative,
							double * prev_x, double * curr_x, VectorD * debug )
{
	*prev_x = *curr_x;
	*curr_x -= evaluate( fnc, *curr_x ) / evaluate( derivative, *curr_x );
	if(debug)
		appendD( debug, *curr_x );
}

double root_newt( const Par_Func * fnc, const Par_Func * derivative, 
					double start, double xtol, double ftol, VectorD * debug )
{
	set_default_tolerance_if_requested( &xtol, &ftol );
	
	if(debug)
	{
		init_debug(debug);
		appendD( debug, start );
	}
	
	double prev_x = 0;
	double curr_x = start;
	bool found_root = false;
	
	for( int i = 0; i < MAX_ITERATION && !found_root; i++ )
	{
		newt_next_step( fnc, derivative, &prev_x, &curr_x, debug );
		found_root = fabs( prev_x - curr_x ) <= xtol 
					&& fabs( evaluate( fnc, curr_x ) ) <= ftol;
	}
	
	if( !found_root )
		raiseErr( "Loop did not stop after %d iterations", MAX_ITERATION );
	
	return curr_x;
}

static inline void sec_next_step( const Par_Func * fnc, double * x_prev , 
	double * x_curr1, double * x_curr2, double * fx_curr1, double * fx_curr2,
			VectorD * debug )
{
	*x_prev   = *x_curr1;
	*x_curr1  = *x_curr2;
	//eprint( "%lf %lf %lf %lf %lf\n", x_prev, x_curr1, x_curr2, fx_curr1, fx_curr2 );
	*x_curr2  = *x_curr1 - *fx_curr2 * ( *x_curr1 - *x_prev ) 
												/ ( *fx_curr2 - *fx_curr1 );
	*fx_curr1 = *fx_curr2;
	*fx_curr2 = evaluate( fnc, *x_curr2 );
	
	
	if(debug)
		appendD( debug, *x_curr2 );
}

double root_sec( const Par_Func * fnc, double x1, double x2, 
					 double xtol, double ftol, VectorD * debug, int * exit_status )
{
	set_default_tolerance_if_requested( &xtol, &ftol );
	
	if(debug)
	{
		init_debug(debug);
		appendD( debug, x1 );
		appendD( debug, x2 );
	}
	
	double x_prev = 0;
	double x_curr1 = x1;
	double x_curr2 = x2;
	double fx_curr1 = evaluate( fnc, x_curr1 );
	double fx_curr2 = evaluate( fnc, x_curr2 );
	bool found_root = false;
	
	//eprint( "%lf %lf %lf %lf %lf\n", x_prev, x_curr1, x_curr2, fx_curr1, fx_curr2 );
	
	for( int i = 0; i < MAX_ITERATION && !found_root; i++ )
	{
		sec_next_step( fnc,  &x_prev , &x_curr1,  &x_curr2,  &fx_curr1,  
		&fx_curr2, debug );
		
		found_root = fabs( x_curr1 - x_curr2 ) <= xtol 
					&& fabs( fx_curr2 ) <= ftol;
	}
	
	if( !found_root )
	{
		if( exit_status )
			*exit_status = -1;
		else
			raiseErr( "Loop did not stop after %d iterations", MAX_ITERATION );
	}
	
	if(exit_status)
		*exit_status = 0;
	
	return x_curr2;
}

/*
double root_inv_qdr( const Par_Func * fnc, double x1, double x2, double x3,
			double xtol, double ftol, VectorD * debug)
{
	set_default_tolerance_if_requested( &xtol, &ftol );
	if(debug)
	{
		init_debug(debug);
		appendD( debug, x1 );
		appendD( debug, x2 );
	}
	
	
	
	
	
	
}

*/



