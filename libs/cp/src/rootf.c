#include "rootf.h"
#include "float.h"
#include <math.h>

#define DEFAULT_XTOLERANCE DBL_EPSILON * 100
#define DEFAULT_FTOLERANCE DBL_EPSILON * 10'000
#define MAX_ITERATION 100'000

static inline int sign( double x )
{
	if( x > 0 )
		return 1;
	if( x < 0 )
		return -1;
	return 0;
}

static inline bool bis_next_step(double *curr, int *sign1, int *sign_curr,
	double *x1, double *x2, ParamFuncPtr fnc, void *params)
{
	*sign_curr = sign( fnc(*curr,params) );
	if( *sign1 == *sign_curr )
	{
		*x1 = *curr;
		*sign1 = *sign_curr;
	} else if( *sign_curr == 0 )
		return true;
	else
		*x2 = *curr;
		
	return false;
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
		curr = x1 + x2 / 2;
		if(bis_next_step(&curr,&sign1,&sign_curr,&x1,&x2,fnc,params))
			return curr;
		if(debug)
			arr_append( debug, curr );
	}
	return curr;
}

static inline void set_default_tolerance_if_requested( double * xtol, double * ftol)
{
	if( *xtol == DEF_TOL )
		*xtol = DEFAULT_XTOLERANCE;
		
	if( *ftol == DEF_TOL )
		*ftol = DEFAULT_FTOLERANCE;
}

static inline void newt_next_step( ParamFuncPtr f, void *params_f,
								   ParamFuncPtr derivative, void *params_der,
							double * prev_x, double * curr_x, Array *debug )
{
//	static int count = 0;
//	count++;
//	if( (count%100) == 0 )
//		eprint("%d\n", count);
	
	*prev_x = *curr_x;
	*curr_x -= f(*curr_x,params_f) / derivative(*curr_x,params_der);
	if(debug)
		arr_append(debug, *curr_x);
		
}

double root_newt( ParamFuncPtr f, void *params_f, 
				  ParamFuncPtr derivative, void *params_der,
				  double start, double xtol, double ftol, Array *debug )
{
	set_default_tolerance_if_requested( &xtol, &ftol );
	
	if(debug)
	{
		arr_free(debug);
		arr_append( debug, start );
	}
	
	double prev_x = 0;
	double curr_x = start;
	bool found_root = false;
	
	int i = 0;
	
	for( ; i < MAX_ITERATION && !found_root; i++ )
	{
		newt_next_step(f,params_f,derivative,params_der, &prev_x, &curr_x, debug );
		found_root = fabs( prev_x - curr_x ) <= xtol 
					&& fabs( f(curr_x,params_f) ) <= ftol;
	}
	
//	printf("%d\n", i);
	
	if( !found_root )
	{
		eprint( "Loop did not stop after %d iterations", MAX_ITERATION );
		curr_x = NAN;
	}
	
	return curr_x;
}

double root_newt_bracket( ParamFuncPtr f, void *params_f, 
				  ParamFuncPtr derivative, void *params_der,
				  double x1,double x2, double xtol, double ftol, Array debug[] )
{
	set_default_tolerance_if_requested( &xtol, &ftol );
	
	if(debug)
	{
		arr_free(debug);
		arr_free(debug+1);
		arr_append( debug, x1);
		arr_append( debug, x2);
	}
	
	bool found_root = false;
	
	double curr = 0;
	double prev_x1=x1;
	int sign1 = sign( f(x1,params_f) );
	int sign_curr = 0;
	
	int i = 0;
	
	for( ; i < MAX_ITERATION && !found_root; i++ )
	{
		prev_x1=x1;
		newt_next_step(f,params_f,derivative,params_der, &prev_x1, &curr,NULL);
		if( x1<curr && curr<x2)
		{
			bis_next_step(&curr,&sign1,&sign_curr,&x1,&x2,f,params_f);
		}else{
			curr = ( x1 + x2 ) / 2;
			bis_next_step(&curr,&sign1,&sign_curr,&x1,&x2,f,params_f);
		}
		
		found_root = fabs( x1 - x2 ) <= xtol 
					&& fabs( f(curr,params_f) ) <= ftol;
					
		arr_append(debug,x1);
		arr_append(debug,x2);
	}
	
//	printf("%d\n", i);
	
	if( !found_root )
	{
		eprint( "Loop did not stop after %d iterations", MAX_ITERATION );
		curr = NAN;
	}
	
	return curr;
}
