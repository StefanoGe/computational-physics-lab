#include "rootf.h"
#include "float.h"
#include <math.h>

#define DEFAULT_XTOLERANCE DBL_EPSILON * 100
#define DEFAULT_FTOLERANCE DBL_EPSILON * 10'000
#define MAX_ITERATION 1000

static inline int sign( double x )
{
	if( x > 0 )
		return 1;
	if( x < 0 )
		return -1;
	return 0;
}

static inline bool bis_next_step(double *curr, int *sign1, int *sign_curr,
	double *x1, double *x2, ParamFuncPtr fnc, void *params, double *fcurr)
{
	*fcurr = fnc(*curr, params);
	*sign_curr = sign( *fcurr );
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
	int sign1 = sign( fnc(x1,params) );
	int sign_curr = 0;
	double curr = 0;
	double fcurr = 0;
	
	if(tol == DEF_TOL)
		tol = DBL_EPSILON * 10;
	tol += DBL_EPSILON * MAX( fabs(x1), fabs(x2) );
	
	int i=0;
	
	for(; i<MAX_ITERATION && fabs( x1 - x2 ) > tol; i++)
	{
		curr = (x1 + x2) / 2;
		if( bis_next_step(&curr,&sign1,&sign_curr,&x1,&x2,fnc,params, &fcurr) )
		{
			if(debug)
				arr_append( debug, curr );
			return curr;
		}
		if(debug)
			arr_append( debug, curr );
	}
	
	if(i >= MAX_ITERATION)
		raiseErr("loop did not terminate after %d iterations", MAX_ITERATION);
	
	return curr;
}

static inline void set_default_tolerance_if_requested( double * xtol, double * ftol)
{
	if( *xtol == DEF_TOL )
		*xtol = DEFAULT_XTOLERANCE;
		
	if( *ftol == DEF_TOL )
		*ftol = DEFAULT_FTOLERANCE;
}

static inline void newt_next_step( ParamFuncPtr derivative, void *params_der,
							double *prev_x, double *curr_x, double *fcurr,
							Array *debug )
{
//	static int count = 0;
//	count++;
//	if( (count%100) == 0 )
//		eprint("%d\n", count);
	
	*prev_x = *curr_x;
	*curr_x -= *fcurr / derivative(*prev_x, params_der);
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
		arr_append(debug, start);
	}
	
	double prev_x = 0;
	double curr_x = start;
	double fcurr = f(curr_x, params_f);
	bool found_root = false;
	
	int i = 0;
	
	for( ; i < MAX_ITERATION && !found_root; i++ )
	{
		newt_next_step(derivative, params_der, &prev_x, &curr_x, 
			&fcurr, debug );
		fcurr = f(curr_x, params_f);
		found_root = fabs( prev_x - curr_x ) <= xtol 
					&& fabs( fcurr ) <= ftol;
	}
	
	if( !found_root )
	{
		eprint( "Loop did not stop after %d iterations", MAX_ITERATION );
		curr_x = NAN;
	}
	
	return curr_x;
}

double root_newt_bracket( ParamFuncPtr f, void *params_f, 
				  ParamFuncPtr derivative, void *params_der,
				  double x1, double x2, double xtol, double ftol, Array *debug)
{
	set_default_tolerance_if_requested( &xtol, &ftol );
	
	bool found_root = false;
	
	double curr = (x2+x1)/2;
	double fcurr = f(curr, params_f);
	double prev=0;
	int sign1 = sign( f(x1,params_f) );
	int sign_curr = 0;
	
	if(debug)
	{
		arr_free(debug);
		arr_append(debug, curr);
	}
	
	int i = 0;
	
	for( ; i<MAX_ITERATION && !found_root; i++ )
	{
		newt_next_step(derivative, params_der, &prev, &curr,
			&fcurr, NULL);
			
		if( x1<curr && curr<x2)
		{
			if(bis_next_step(&curr, &sign1, &sign_curr, &x1, &x2, f, params_f,
				&fcurr))
			{
				if(debug)
					arr_append(debug, curr);
				return curr;
			}
		}else{
			curr = ( x1 + x2 ) / 2;
			if(bis_next_step(&curr, &sign1, &sign_curr, &x1, &x2, f, params_f,
				&fcurr))
			{
				if(debug)
					arr_append(debug, curr);
				return curr;
			}
		}
					
		found_root = (fabs(curr - prev) <= xtol || fabs(x2 - x1) <= xtol)
			&& fabs(fcurr) <= ftol;

		
		if(debug)
			arr_append(debug,curr);
	}
	
//	printf("%d\n", i);
	
	if( !found_root )
	{
		eprint( "Loop did not stop after %d iterations", MAX_ITERATION );
		curr = NAN;
	}
	
	return curr;
}

static inline void sec_next_step( const ParamFuncPtr f, void *params,
	double *x_prev,
	double * x_curr1, double * x_curr2, double * fx_curr1, double * fx_curr2,
			Array *debug )
{
	*x_prev   = *x_curr1;
	*x_curr1  = *x_curr2;
	//eprint( "%lf %lf %lf %lf %lf\n", x_prev, x_curr1, x_curr2, fx_curr1, fx_curr2 );
	*x_curr2  = *x_curr1 - *fx_curr2 * ( *x_curr1 - *x_prev ) 
												/ ( *fx_curr2 - *fx_curr1 );
	*fx_curr1 = *fx_curr2;
	*fx_curr2 = f(*x_curr2,params);
	
	
	if(debug)
		arr_append( debug, *x_curr2 );
}

double root_sec( const ParamFuncPtr f, void *params,
					double x1, double x2, 
					double xtol, double ftol, Array *debug, int *exit_status )
{
	set_default_tolerance_if_requested( &xtol, &ftol );
	
	if(debug)
	{
		arr_append( debug, x1 );
		arr_append( debug, x2 );
	}
	
	double x_prev = 0;
	double x_curr1 = x1;
	double x_curr2 = x2;
	double fx_curr1 = f(x_curr1,params);
	double fx_curr2 = f(x_curr2,params);
	bool found_root = false;
	
	//eprint( "%lf %lf %lf %lf %lf\n", x_prev, x_curr1, x_curr2, fx_curr1, fx_curr2 );
	
	for( int i = 0; i < MAX_ITERATION && !found_root; i++ )
	{
		sec_next_step( f, params, &x_prev , &x_curr1,  &x_curr2,  &fx_curr1,  
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



void inv_quad_next_step(ParamFuncPtr f, void *params,
				double *xs, double *fxs, int *older_x_index,
				int *second_x_index, int *newer_x_index)
{	
	int temp=*older_x_index;
	*older_x_index=*second_x_index;
	*second_x_index=*newer_x_index;
	*newer_x_index=temp;
	
	
	for(int i=0; i<3; i++)
	{
		if(fxs[i]==0)
		{
			xs[*newer_x_index]=xs[i];
			fxs[*newer_x_index]=0;
			return;
		}
	}
	
	double weights[3]={1,1,1};
	
	for(int i=0; i<3; i++)
		for(int j=0; j<3; j++)
			if(i!=j)
			{
				weights[i]/=(fxs[i]-fxs[j]);
//				eprint("fxs[%d]-fxs[%d] = %lf",i,j,fxs[i]-fxs[j]);
//				eprint("weights[%d] = %lf", i,weights[i]);
			}
/*
	eprint("%lf %lf %lf\n",weights[*older_x_index],
			weights[*second_x_index],weights[*newer_x_index]);
*/
	double den = 0;
	double num = 0;
	
	for(int i=0; i<3; i++)
	{
		den+=weights[i]/fxs[i];
		num+=weights[i]*xs[i]/fxs[i];
	}
	
//	eprint("num %lf den %lf\n",num,den);
	
	xs[*newer_x_index]=num/den;
	fxs[*newer_x_index]=f(xs[*newer_x_index],params);
	
//	mypause();
}

double root_inv_quad_int(const ParamFuncPtr f, void *params, double x1,
				double x2, double x3, double xtol, double ftol, Array *debug,
				int *exit_status)
{
	set_default_tolerance_if_requested( &xtol, &ftol );
	
	if(debug)
	{
		arr_append( debug, x1 );
		arr_append( debug, x2 );
		arr_append( debug, x3);
	}
	
	double xs[3]={x1,x2,x3};
	double fxs[3]={f(x1,params),f(x2,params),f(x3,params)};
	int older_x_index=0;
	int second_x_index=1;
	int newer_x_index=2;
	bool found_root=0;
	
	for( int i = 0; i < MAX_ITERATION && !found_root; i++ )
	{
		inv_quad_next_step(f,params,xs,fxs,&older_x_index,
			&second_x_index,&newer_x_index);
		
		arr_append(debug,xs[newer_x_index]);
				
		found_root = fabs( xs[second_x_index] - 
						xs[newer_x_index] ) <= xtol 
					&& fabs( fxs[newer_x_index] ) <= ftol;
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
	
	return xs[newer_x_index];
}



